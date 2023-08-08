import matplotlib.pyplot as plt
import csv
import numpy as np

noeud_en_ms = 0.51444468

def calcul_puissance(fichier_vent, fichier_polaires, vnavire, affichage, nomsolution="", retour=False):
    """
    Cette fonction renvoie une puissance en kW pour une route une solution vélique et une vitesse de avire.

    Cette foction renvoie pour un fichier de vent donné (statistiques de vent pour un trajet), un fichier de polaire
    pour une solution vélique, et une vitesse navire (en noeuds); une puissance moyenne sur le trajet en kW.
    Elle permet également d'afficher un graphique des puissances selon les vitesses et l'orientation du vent pour
    connaître le potentiel de la solution pour ce vent.

    :param fichier_vent: type(str) fichier .csv avec séparateur ";" avec les angles sur la première ligne en degrés, et
    les vitesses en m/s sur la première colonne, pour chacun de ceux-ci il y a une probabilité associée.

    :param fichier_polaires: type(str) fichier .csv avec séparateur ";" avec en première ligne les vitesses en noeuds,
    ensuite par ligne, on a une première colonne l'effort en kN, et la seconde l'angle associé pour la vitesse donnée en
    degrés, ce pour chacune des vitesses, il y aura donc 2 fois plus de valeurs pour les lignes suivant la ligne des
    vitesses.

    :param vnavire: type(float) vitesse du navire en nœuds.

    :param affichage type(bool) True ou False pour dire si on affiche ou non le graphique d'une part pour chacune des
    vitesses et d'autre part pour toutes les vitesses combinées.

    :param nomsolution type(str) Nom de la solution étudiée

    :return : La puissance procurée par les voiles pour la route du navire et pour la vitesse donnée
    (important donc d'avoir une vitesse en nœuds)

    """
    if retour:
        v_lecture=0
    else:
        v_lecture=vnavire
    stats_dict = lecture_vent(fichier_vent, v_lecture)  # On stocke les valeursdes statistiques dans un dictionnaire
    if retour:
        stats_retour=creer_stats_retour(stats_dict)
        stats_retour_apparent=conversion_donnée(stats_retour,vnavire)
        stats_dict=stats_retour_apparent
    # {vitesse1:{angle1:proba1, angle2:proba2...}, vitesse2{...},...}
    les_vitesses_stats = list(stats_dict.keys())  # On récupère toutes les vitesses en m/s des stats
    les_angles_stats = list(stats_dict[les_vitesses_stats[0]].keys())  # On récupère les angles des stats
    polaires = lecture_pol(fichier_polaires)  # On stocke les valeurs des polaires dans un dictionnaire
    # {vitesse1:{angle1:effort1, angle2:effort2..}, vitesse2{...},...}
    les_vitesses_polaires = list(polaires.keys())  # On récupère les vitesses des polaires
    proba_effort = 0  # l'effort moyen que va fournir la solution vélique moyenne pondérée par la proba du vent pour
    # une vitesse et un angle
    somme_proba = 0  # Somme de la proba de tous les angles et toutes les vitesses considérées pour voir pendant
    # combien de temps la voile sera utilisable
    affichage2 = True
    bool_v_utilisable = True  # Permet de tester si dans nos boucles on dépasse la vitesse max des polaires
    résultat_par_incidence = {}  # Stocke les efforts sommés par incidence
    for angle in les_angles_stats:
        résultat_par_incidence[angle] = 0  # initialisation
    #########################Calcul de l'effort moyen pour le vent de la route donnée#################################
    for vitesse in les_vitesses_stats:
        # On parcourt les vitesses de vents dont on connait la répartition
        # On cherche à déterminer la vitesse la plus proche dans les polaires en dessous de la vitesse du vent qui
        # permet d'être conservateur et de considéré un effort plus faible que la réalité et servira de vitesse
        # considérée pour le calcul
        borne_sup = vitesse
        vitesse_cal = 0
        les_angles_pol = np.array([])
        for vitesse_pol in les_vitesses_polaires:
            # Boucle qui permet de renvoyer avec la variable vitesse_cal la plus grande valeur inférieure à la vitesse
            # de vent
            if vitesse_pol <= borne_sup and vitesse_pol >= vitesse_cal:
                vitesse_cal = vitesse_pol
            if vitesse > max(les_vitesses_polaires):
                # On ne calculera pas si on dépasse la vitesse des polaires puisqu'on dira qu'à ce moment la solution
                # n'est plus utilisable
                bool_v_utilisable = False
            if vitesse <= max(les_vitesses_polaires):
                bool_v_utilisable = True
        if vitesse_cal != 0 and bool_v_utilisable:
            # On ne peut considérer une vitesse nulle et une vitesse supérieure au max des vitesses des polaires qui
            # rendrait la solution inutilisable
            polaire_vitesse_cal = polaires[vitesse_cal]
            les_angles_pol = list(polaire_vitesse_cal.keys())  # On récupère les angles des polaires pour la vitesse
            # donnée
            save_les_angles_pol = np.copy(les_angles_pol)  # On sauvegarde les angles des polaires pour pouvoir en
            # générer la symétrie des valeurs
            les_efforts = []  # Récupération des efforts associée à la vitesse_cal
            for angle in les_angles_pol:
                les_efforts.append(polaire_vitesse_cal[angle])
            # Symétrisation des efforts pour les polaires
            for i in range(-1, -len(save_les_angles_pol) - 1, -1):
                if save_les_angles_pol[i] != 180:
                    les_angles_pol.append(360 - save_les_angles_pol[i])
                    les_efforts.append(polaire_vitesse_cal[save_les_angles_pol[i]])
                else:
                    pass
            les_efforts = np.array(les_efforts)
            # On a donc nos efforts pour vitesse_cal et rangé de 0 à 359.9° dans les_efforts
            interp_efforts = np.interp(les_angles_stats, np.array(les_angles_pol), les_efforts)  # On veut connaitre la
            # valeur de nos efforts pour les angles de nos vents statistiques on interpole donc pour x=les_angles_pol
            # et f(x)=les_efforts
            effort_dict = {}
            # On range tous les résultats dans un dico effort_dict
            for i in range(len(les_angles_stats)):
                effort_dict[les_angles_stats[i]] = interp_efforts[i]
            # Pour les angles stats on calcule l'effort associé à partir de sa proba comme une moyenne pondérée
            for angle in les_angles_stats:
                effort = effort_dict[angle]
                proba_effort += effort * stats_dict[vitesse][angle]  # Pour chacune des vitesses on ajoute
                # l'effort pour la direction
                somme_proba += stats_dict[vitesse][angle]  # Ajout de la proba considérée
                résultat_par_incidence[angle] += effort * stats_dict[vitesse][angle] * vnavire * noeud_en_ms
                # On ajoute pour chacune des vitesses l'effort proba pour l'angle au dico résultat par incidence
            if affichage:
                # Si affichage == True on pourra afficher un graphique avec ces valeurs pour chacune des vitesses suivant
                # l'orientation du vent
                les_angle = list(effort_dict.keys())
                les_r = [effort_dict[angle] * stats_dict[vitesse][angle] * vnavire * noeud_en_ms for angle in les_angle]
                les_angle = les_angle[::-1]  # inversion de la liste pour passer en sens horaire, mais pas des valeurs
                les_r = np.array(les_r)
                les_angle = np.array(les_angle) * np.pi / 180
                plt.polar(les_angle, les_r, label=str(vitesse))
                etiquettes_angles = np.linspace(45, 315, 7)
                etiquettes_angles = np.append(etiquettes_angles, 0)
                etiquettes_angles = etiquettes_angles[::-1]
                plt.gca().set_theta_offset(np.pi / 2)
                plt.gca().set_xticklabels(etiquettes_angles)
                plt.subplots_adjust(top=0.8)
    if affichage:
        plt.legend(loc='upper right', bbox_to_anchor=(1.37, 1.0))
        plt.title(
            "Puissance (en kW) pondérée par la probabilité du vent suivant\n l'orientation du vent pour différentes "
            "vitesses du vent (en m/s)")
        plt.show()
    if affichage:
        # Si affichage==True on pourra afficher un graphique de l'orientation des efforts suivant le vent
        les_angles = list(résultat_par_incidence.keys())
        les_dists = []
        for x in les_angles:
            les_dists.append(résultat_par_incidence[x])
        les_angles = les_angles[::-1]  # inversion de la liste pour passer en sens horaire mais pas des valeurs qui
        # elles restent dans le bon sens
        les_angles = np.array(les_angles)
        les_angles = np.pi * les_angles / 180
        les_dists = np.array(les_dists)
        les_r = les_dists
        plt.figure(figsize=(8, 8))
        plt.polar(les_angles, les_r, label=str(vnavire))
        etiquettes_angles = np.linspace(45, 315, 7)
        etiquettes_angles = np.append(etiquettes_angles, 0)
        etiquettes_angles = etiquettes_angles[::-1]
        plt.gca().set_theta_offset(np.pi / 2)
        plt.gca().set_xticklabels(etiquettes_angles)
        plt.subplots_adjust(top=0.8)
        plt.legend()
        plt.title("Puissance suivant l'incidence pondérée par la probabilité de l'incidence pour " + nomsolution)
        plt.show()
    # avec un test sur excel on a pas oublié d'angle ni de vitesse
    return proba_effort * vnavire * noeud_en_ms, résultat_par_incidence


def comparaison_techno(vitesse, list_techno, stats_vent, affichage, retour=False):
    """
    Fonction qui permet de comparer les solutions de voiles entre elles et renvoie la plus optimisée.

    Cette fonction prend une vitesse de navire donnée ainsi qu'une liste de technologie et des informations de
    statistiques de vent pour une route, elle va faire appel à la fonction calcul_puissance pour chacune des solutions,
    et va renvoyer celle pour laquelle la solution est la plus optimisée.

    :param vitesse: type(float), vitesse du navire en noeuds, doit être la même que celle du fichier de polaire de
    solution vélique
    :param list_techno: type (list) Liste de chaine de caractères correspondant aux solutions à étudier
    "" "polaires/" + techno + "_" + str(vitesse) + ".csv" "" chemin utilisé pour trouver le fichier associé à la techno
    :param stats_vent: type(str) fichier .csv avec séparateur ";" avec les angles sur la première ligne en degrés, et
    les vitesses en m/s sur la première colonne, pour chacun de ceux-ci il y a une probabilité associée
    :param affichage type(bool) Permet d'afficher un graphique avec le résultat pour chacune des solutions avec les
    vents combinés (puissance dépendant seulement de l'orientation et pas de la vitesse)
    :return: Le nom de la technologie la plus optimisée pour la route donnée, ainsi que la puissance que celle-ci fait
    économiser
    """
    res = {}  # Stocke les résultats clé puissance et techno en valeur
    for techno in list_techno:
        nom_fichier = "polaires/" + techno + "_" + str(11) + ".csv"
        res_techno, dico = calcul_puissance(stats_vent, nom_fichier, vitesse, False, techno, retour)
        # On récupère le résultat pour la techno puissance moyenne pondérée et dictionnaire associé pour les différentes
        # incidences
        if affichage:
            # Si affichage==True on affiche pour chacune des technos une courbe polaire pour la puissance suivant
            # l'incidence
            les_angles = list(dico.keys())
            les_val = []
            for angle in les_angles:
                les_val.append(dico[angle])
            les_angles = les_angles[::-1]  # inversion de la liste pour passer en sens horaire, on inverse pas
            # le sens des valeurs
            les_angles = np.radians(les_angles)
            plt.polar(les_angles, les_val, label=techno)
            etiquettes_angles = np.linspace(45, 315, 7)
            etiquettes_angles = np.append(etiquettes_angles, 0)
            etiquettes_angles = etiquettes_angles[::-1]
            plt.gca().set_theta_offset(np.pi / 2)
            plt.gca().set_xticklabels(etiquettes_angles)
            plt.subplots_adjust(top=0.8)
        res[res_techno] = techno
    if affichage:
        plt.legend(loc='upper right', bbox_to_anchor=(1.37, 1.0))
        plt.title(
            "Comparaison des puissances fournies (en kW) selon l'orientation\n du vent pondérée par la probabilité "
            "de l'orientation de celui-ci")
        plt.show()
    vals = list(res.keys())
    opti = max(vals)
    techno_opti = res[opti]
    print(res,stats_vent)
    return (techno_opti, opti, res)


def comparaison_route(vitesse: object, list_techno: object, liste_route: object) -> object:
    """Cette fonction compare les différentes routes pour les différentes technos, et renvoie la techno la plus
    optimisée pour la route la plus optimisée

    :param vitesse: type(float), vitesse du navire en noeuds, doit être la même que celle du fichier de polaire de
    solution vélique
    :param list_techno: type (list) Liste de chaine de caractères correspondant aux solutions à étudier
    "" "polaires/" + techno + "_" + str(vitesse) + ".csv" "" chemin utilisé pour trouver le fichier associé à la techno
    :param liste_route: type (list) Liste des chemins des différentes routes

    :return: route optimisée, techno optimale pour la route, puissance moyenne pour cette route et cette techno
    """
    res = {} #Stockage des valeurs pour la route clé route et valeur le max pour la route considérée
    for route in liste_route:
        dictionnaire = comparaison_techno_aller_retour(vitesse, list_techno, route)
        les_puissances=list(dictionnaire.keys())
        max_puissance= max(les_puissances)
        res[route] = (dictionnaire[max_puissance],max_puissance )
    #on stocke la meilleure solution pour chacune des routes
    vals = list(res.keys())
    maximum = 0
    for valeur in vals:
        résultat = res[valeur]
        if résultat[1] > maximum:
            # Récupération du meilleur résultat pour les différentes routes
            maximum = résultat[1]
            techno_opti = résultat[0]
            route_opti = valeur
        else:
            pass
    print(techno_opti + " me fait économiser " + str(maximum) + " par la route " + route_opti)
    return route_opti, techno_opti, maximum


def affichage_comparaison_route(liste_route,vnavire=0, retour=False):
    """Cette fonction permet d'afficher les vents de différentes routes.

    :param liste_route: type (list) Liste des chemins des différentes routes
    :return: Affiche un graphique avec les stats des vents pour différentes routes en combinant la proba des
    différentes vitesses
    """
    for i in range(len(liste_route)):
        route = liste_route[i]
        if i < len(liste_route) - 1:
            end = False
        else:
            end = True
        if retour:
            # On affiche un graphique qui compare les différentes routes si on choisit d'afficher ou pour une seule
            # si telle est son utilisation
            dico = lecture_vent(route, 0)
            dico_retour=creer_stats_retour(dico)
            if vnavire!=0:
                stats_retour = conversion_donnée(dico_retour, vnavire)
            else:
                stats_retour=dico_retour
        else:
            dico = lecture_vent(route, vnavire)
            stats_retour = dico
        probas={}
        les_vitesses=list(stats_retour.keys())
        for angle in stats_retour[les_vitesses[0]]:
            probas[angle]=0
        for vitesse in stats_retour:
            for angle in stats_retour[vitesse]:
                probas[angle]+=stats_retour[vitesse][angle]
        les_r = list(probas.values())
        les_angles = list(probas.keys())
        les_angles = les_angles[
                     ::-1]  # inversion de la liste pour passer en sens horaire, les_r restent eux dans le
        # bon ordre ca a été vérifié par affichage
        plt.polar(np.radians(les_angles), les_r, label=route)
        etiquettes_angles = np.linspace(45, 315, 7)
        etiquettes_angles = np.append(etiquettes_angles, 0)
        etiquettes_angles = etiquettes_angles[::-1]
        plt.gca().set_theta_offset(np.pi / 2)
        plt.gca().set_xticklabels(etiquettes_angles)
        plt.subplots_adjust(top=0.8)
        plt.title("Probabilité d'occurrence des vents suivant l'orientation")

    plt.legend(loc='upper right', bbox_to_anchor=(1.37, 1.0))
    plt.show()
    return


def lecture_pol(fichier):
    """
    Cette fonction permet de lire le fichier des polaires associé à une solution vélique et renvoie un dictionnaire
    {vitesse1:{angle1:effort1, angle2:effort2..}, vitesse2{...},...
    :param fichier: type(str) Correspondant au chemin d'accès au fichier des polaires (vitesse en noeuds)
    :return: dictionnaire {vitesse1:{angle1:effort1, angle2:effort2..}, vitesse2{...},...}
    """
    polaires = open(fichier, "r", encoding="utf-8")  # Ouverture de fichier de polaire
    polaires_lines = csv.reader(polaires)
    polaires = {}  # On stocke les polaires sous forme d'un dictionnnaire {vitesse1:{angle1:effort1,angle2:effort2,
    # ...},vitesse2:{...},...}
    indentation = 0
    les_vitesses_polaires = []  # Pour stocker les vitesses des polaires
    # Attention les vitesses polaires sont en noeuds dans notre fichier on les converties en m/s
    save = []  # On sauvegarde pour plus tard les lignes dans une liste
    probas = {}
    for line in polaires_lines:
        line_formatted_inter = line[0].strip().split(";")  # Mise en forme
        line_formatted = []
        for val in line_formatted_inter:
            try:
                if len(val) == 0:
                    pass
                else:
                    line_formatted.append(float(val))  # On récupère les valeurs dans une liste
            except:
                line_formatted.append(val)  # Si on a des caractères on les garde quand même
        if indentation == 0:
            for val2 in line_formatted:
                if type(val2) == float:  # Pour éviter les chaines de caractère
                    les_vitesses_polaires.append(val2 * noeud_en_ms)  # Stockage des vitesses en ms pour les
                    # comparer à celles des stats
        else:
            save.append(line_formatted)  # Sauvegarde des lignes
        indentation += 1  # Incrémentation de l'indentation
    ###################################################################################################################
    for vitesse in les_vitesses_polaires:
        polaires[vitesse] = {}  # Pour chacune des vitesses on a des polaires {angle1:effort1, angle2:effort2,...}
    # Dans la liste vals on a [[effort_vitesse1_1,angle_vitesse1_1,effort_vitesse2_1,angle_vitesse2_1,...],
    # [effort_vitesse1_2,angle_vitesse1_2,effort_vitesse2_2,angle_vitesse2_2,...],...].
    # Car dans chacune des lignes on a les efforts les angles rangés dans l'ordre des vitesses
    ##########################################Génération du dictionaire de polaires###################################
    for vals in save:
        n = len(vals)  # Nombre total de valeurs dans la ligne ce qui équivaut à 2 fois le nombre de vitesse
        # (car un angle et un effort pour chaque vitesse)
        effort = 0
        angle = 0
        for i in range(n):
            if i % 2 == 0:
                # Quand i est pair on est dans la première colone associée à la vitesse donc la force (en kN)
                effort = vals[i]
            else:
                # Quand i est impaire on est dans la deuxième donc on a la valeur de l'angle associée
                angle = vals[i]
                polaires[les_vitesses_polaires[i // 2]][angle] = effort  # Ajout de l'effort pour la vitesse associée
                # et l'angle associé
                effort = 0
                angle = 0
    return polaires


def lecture_vent(fichier, vnavire=0, affichage=False, end=True):
    """
    Cette fonction permet de lire un fichier vent réel et renvoie les probabilités d'obtenir un vent apparent. :param
    fichier: type(str) chemin d'accès au fichier de vent (vitesse en m/s) :param vnavire: type(float) vitesse du
    navire (en noeuds), par défaut nulle :param affichage: type(bool) permet d'afficher ou non un graphique comparant
    les différentes vitesses de navire :param end: type(bool) Permet lorsque utilisé dans la fonction
    affichage_comparaison_route pour le dernier, on affiche la figure avec toutes les courbes :return: Type(dic) {
    vitesse1:{angle1:proba1, angle2:proba2,...}, vitesse2:{...},...} Et affiche un graphique si affichage==True
    """
    stats = open(fichier, "r", encoding="utf-8")  # ouverture du fichier des stats vent
    stats_lines = csv.reader(stats)
    stats_dict = {}  # On stocke les valeurs dans un dictionnaire {vitesse1:{angle1:proba1, angle2:proba2...},
    # vitesse2{...},...}
    indentation = 0  # utile pour savoir dans quelle ligne on est
    les_angles_stats = []  # liste des angles
    les_vitesses_stats = []  # liste des vitesses
    probas = {}
    #######################################Lecture du fichier de stats#################################################
    for line in stats_lines:
        line_formatted = line[0].strip().split(";")  # mise en forme
        line_formatted1 = []
        for val in line_formatted:
            if len(val) == 0:
                pass  # on ne récupère pas les caractères vides
            else:
                val2 = float(val)
                line_formatted1.append(val2)  # les valeurs dans une liste
        if indentation == 0:
            les_angles_stats = line_formatted1  # dans la première ligne se trouve les angles
        if len(les_angles_stats) != 0:
            stats_dict[line_formatted1[0]] = {}  # On prend un dictionnaire pour cette vitesse car la vitesse est
            # stockée dans la première colonne
            for i in range(len(line_formatted1)):
                if i != 0:
                    stats_dict[line_formatted1[0]][les_angles_stats[i - 1]] = line_formatted1[i]  # On ajoute les
                    # valeurs pour chacun des angles associés
                else:
                    les_vitesses_stats.append(line_formatted1[i])
        indentation += 1  # incrémentation d'indentation pour dire que l'on passe à la ligne suivante
    for angle in les_angles_stats:
        probas[angle]=0
    if vnavire != 0:
        dic_vnavire = conversion_donnée(stats_dict, vnavire)
        for vitesse in dic_vnavire:
            for angle in dic_vnavire[vitesse]:
                probas[angle] += dic_vnavire[vitesse][angle]
                # Stocke la somme des probas pour les différentes vitesses pour l'orientation d'un vent
    elif vnavire==0:
        for vitesse in stats_dict:
            for angle in stats_dict[vitesse]:
                probas[angle] += stats_dict[vitesse][angle]
                # Stocke la somme des probas pour les différentes vitesses pour l'orientation d'un vent
        #Permet d'avoir le vent apparent
    for angle in les_angles_stats:
        probas[angle] = 0
    if vnavire!=0:
        return dic_vnavire
    elif vnavire==0:
        return stats_dict


def économie_puissance(puissance_sans_voile, puissance_avec_voile, stats_vent, vnavire, vvent, affichage=False):
    pol_voil = lecture_pol(puissance_avec_voile)
    vvent = vvent * noeud_en_ms
    vnavire = vnavire * noeud_en_ms
    pol_sans_voil = lecture_pol(puissance_sans_voile)
    stats_vent = lecture_vent(stats_vent, vnavire)
    puissance_totale = 0
    for vitesse in pol_sans_voil.keys():
        if vitesse == vnavire:
            val_without_s = {}
            val_s = {}
            angle_without_s = list(pol_sans_voil[vitesse].keys())
            angle_s = list(pol_voil[vitesse].keys())
            windspeed_calc = 0
            for windspeed_stat in stats_vent:
                if windspeed_stat >= windspeed_calc and windspeed_stat < vvent:
                    windspeed_calc = windspeed_stat
            stats = stats_vent[windspeed_calc]
            les_angles = list(stats.keys())
            res = {}
            for angle in angle_without_s:
                if angle == 0 or angle == 180:
                    val_without_s[angle] = pol_sans_voil[vitesse][angle]
                else:
                    val_without_s[360 - angle] = pol_sans_voil[vitesse][angle]
                    val_without_s[angle] = pol_sans_voil[vitesse][angle]
            for angle in angle_s:
                if angle == 180:
                    val_s[angle] = pol_voil[vitesse][angle]
                else:
                    val_s[360 - angle] = pol_voil[vitesse][angle]
                    val_s[angle] = pol_voil[vitesse][angle]
            sort_angle_without_s = np.array(list(val_without_s.keys()))
            sort_angle_without_s = np.sort(sort_angle_without_s)
            sort_angle_s = np.array(list(val_s.keys()))
            sort_angle_s = np.sort(sort_angle_s)
            sort_val_without_s = [val_without_s[angle] for angle in sort_angle_without_s]
            sort_val_s = [val_s[angle] for angle in sort_angle_s]
            interp_without_s = np.interp(les_angles, sort_angle_without_s, sort_val_without_s)
            interp_s = np.interp(les_angles, sort_angle_s, sort_val_s)
            les_probas = list(stats.values())
            somme = 0
            test = []
            for proba in les_probas:
                somme += proba
            for i in range(len(les_angles)):
                angle = les_angles[i]
                res[angle] = (interp_without_s[i] - interp_s[i]) * stats[angle] / somme
                test.append((interp_without_s[i] - interp_s[i]) * stats[angle] / somme)
                puissance_totale += (interp_without_s[i] - interp_s[i]) * stats[angle] / somme
            if affichage:
                les_angles = les_angles[::-1]  # inversion de la liste pour passer en sens horaire
                plt.polar(np.radians(les_angles), test, label=str(vvent))
                plt.gca().set_theta_offset(np.pi / 2)
                etiquettes_angles = np.linspace(45, 315, 7)
                etiquettes_angles = np.append(etiquettes_angles, 0)
                etiquettes_angles = etiquettes_angles[::-1]
                plt.gca().set_xticklabels(etiquettes_angles)
                plt.subplots_adjust(top=0.8)
                plt.legend(loc='upper right', bbox_to_anchor=(1.37, 1.0))
                plt.title("Puissance (en kW) pondérée par la probabilité du vent suivant\n l'orientation du vent pour "
                          "différentes vitesses du vent (en m/s)")
                plt.show()
    print(puissance_totale)
    return res


def arrondir_a_5(val):
    """
    Cette fonction permet d'arrondir à une valeur à 5 près
    :param val: type(float)
    :return: type(int)
    """
    quotient = val / 5
    arrondi = round(quotient)
    resultat = arrondi * 5
    return resultat


def conversion_donnée(dico, vnavire):
    """
    Fonction qui convertit les données de vent réel en un vent apparent pour une vitesse navire donnée en noeuds.

    :param dico: type(dictionnaire) Dictionnaire de statistique des vents {vitesse1:{angle1:proba1, angle2:proba2,...},
    vitesse2:{...},...} vitesses en m/s
    :param vnavire: type(float) vitesse en noeuds
    :return: type dictionnaire
    """
    dico_act = {}   # Dictionnaire actualisé pour un vent apparent, organisé de la même manière que pour les
    # statistiques vent du lecture vent, {vitesse1:{angle1:proba1,angle2:proba2,...},vitesse2{...},...}
    vnavire_ms = vnavire * noeud_en_ms
    les_vitesses = list(dico.keys())
    les_angles = list(dico[0].keys())
    for vitesse in les_vitesses:
        dico_act[vitesse] = {}
    somme_rest = 0
    for vitesse_vent in dico:
        for angle_vent in dico[vitesse_vent]:
            angle_réel_vent = angle_vent
            x, y = vitesse_vent * np.cos(angle_réel_vent * np.pi / 180), vitesse_vent * np.sin(
                angle_réel_vent * np.pi / 180)
            vect_v_vent = np.array([x, y])
            vect_v_navire = np.array([vnavire_ms, 0])
            vect_v_réel = vect_v_vent + vect_v_navire
            angle_réel = np.arctan2(vect_v_réel[1], vect_v_réel[0])
            angle_calc=angle_réel
            if angle_réel < 0 and -np.pi <= angle_réel:
                angle_calc =angle_réel + 2 * np.pi
            angle_calc *= 180 / np.pi
            # angle_réel=conversion_angle_repère_bat(angle_réel)
            v_réelle = np.linalg.norm(vect_v_réel)
            v_réelle_ar = round(v_réelle)
            if v_réelle_ar in les_vitesses:
                if angle_calc in dico_act[v_réelle_ar]:
                    dico_act[v_réelle_ar][angle_calc] += dico[vitesse_vent][angle_vent]
                else:
                    dico_act[v_réelle_ar][angle_calc] = dico[vitesse_vent][angle_vent]
            if v_réelle_ar not in les_vitesses:
                somme_rest += dico[vitesse_vent][angle_vent]
    dico_fin = {}
    for vitesse in les_vitesses:
        dico_fin[vitesse] = {}
        for angle in les_angles:
            dico_fin[vitesse][angle] = 0
    somme_proba = 0
    for vitesse in dico_act:
        for angle in dico_act[vitesse]:
            angle_ar = arrondir_a_5(angle)
            if angle_ar == 360:
                angle_ar = 0
            dico_fin[vitesse][angle_ar] += dico_act[vitesse][angle]
            somme_proba += dico_act[vitesse][angle]
    return dico_fin


def affichage_des_polaires(dico):
    for vitesse in dico:
        les_angles = list(dico[vitesse].keys())
        les_angles = np.array(les_angles)
        les_angles = 360 - les_angles
        les_r = list(dico[vitesse].values())
        etiquettes_angles = np.linspace(45, 315, 7)
        etiquettes_angles = np.append(etiquettes_angles, 0)
        etiquettes_angles = etiquettes_angles[::-1]
        plt.polar(np.radians(les_angles), les_r, label=str(round(vitesse, 2)))
        plt.gca().set_theta_offset(np.pi / 2)
        plt.gca().set_xticklabels(etiquettes_angles)
        plt.subplots_adjust(top=0.8)
        plt.title("Polaire de la technologie pour les différentes vitesses (en m/s)")
        plt.legend(loc='upper right', bbox_to_anchor=(1.37, 1.0))
    plt.show()
    return


def creer_stats_retour(dico):
    dico_act={}
    for vitesse in dico:
        dico_act[vitesse]={}
        for angle in dico[vitesse]:
            if angle<180:
                angle_act=180+angle
            elif angle>=180:
                angle_act=angle-180
            dico_act[vitesse][angle]=dico[vitesse][angle_act]
    return dico_act


def comparaison_techno_aller_retour(vitesse: object, list_techno: object, route):
    technos_aller=comparaison_techno(vitesse,list_techno,route,True,False)
    technos_retour=comparaison_techno(vitesse,list_techno,route,True,True)
    dico_aller=technos_aller[2]
    dico_retour=technos_retour[2]
    moy={}
    for val_aller in dico_aller:
        for val_retour in dico_retour:
            if dico_aller[val_aller]==dico_retour[val_retour]:
                val=(val_retour+val_aller)/2
                moy[val]=dico_retour[val_retour]
    return moy


def calc_puissance_pour_vitesse(fichier_stat,fichier_polaire,vitesse_vent,v_navire):
    dic_vent=lecture_vent(fichier_stat,11)
    dic_polaire=lecture_pol(fichier_polaire)
    v_vent_ms=vitesse_vent*noeud_en_ms
    v_stats=list(dic_vent.keys())
    v_polaire=list(dic_polaire.keys())
    v_vent_calc=0
    for v in v_stats:
        if v<=v_vent_ms and v>=v_vent_calc:
            v_vent_calc=v
    v_pol_calc=0
    for v in v_polaire:
        if v<=v_vent_calc*1.1 and v>=v_pol_calc:
            v_pol_calc=v
    proba_effort=0
    s_proba=0
    les_angles_stats=list(dic_vent[v_vent_calc])
    les_angles_pol=list(dic_polaire[v_pol_calc])
    for i in range (-1,-(len(les_angles_pol)-1),-1):
        angle=les_angles_pol[i]
        if angle!=180 and angle!=0:
            dic_polaire[v_pol_calc][360-angle]=dic_polaire[v_pol_calc][angle]
    les_angles_pol=list(dic_polaire[v_pol_calc])
    les_efforts=list(dic_polaire[v_pol_calc].values())
    interp_efforts = np.interp(les_angles_stats, np.array(les_angles_pol), les_efforts)  # On veut connaitre la
    val_interp={}
    les_r={}
    for i in range (len(les_angles_stats)):
        val_interp[les_angles_stats[i]]=interp_efforts[i]
    for angle in dic_vent[v_vent_calc]:
        proba_effort+=dic_vent[v_vent_calc][angle]*val_interp[angle]*v_navire*noeud_en_ms
        s_proba+=dic_vent[v_vent_calc][angle]
        les_r[angle]=dic_vent[v_vent_calc][angle]*val_interp[angle]*v_navire*noeud_en_ms
    etiquettes_angles = np.linspace(45, 315, 7)
    etiquettes_angles = np.append(etiquettes_angles, 0)
    etiquettes_angles = etiquettes_angles[::-1]
    list_r=[]
    s=0
    for angle in dic_vent[v_vent_calc]:
        list_r.append(les_r[angle]/s_proba)
        s+=les_r[angle]/s_proba
    print(s)
    plt.polar(np.radians(les_angles_stats), list_r, label=str(round(v_vent_ms, 2)))
    plt.gca().set_theta_offset(np.pi / 2)
    plt.gca().set_xticklabels(etiquettes_angles)
    plt.subplots_adjust(top=0.8)
    plt.title("Polaire de la technologie pour les différentes vitesses (en m/s)")
    plt.legend(loc='upper right', bbox_to_anchor=(1.37, 1.0))
    plt.show()
    print("proba vent", s_proba)
    return proba_effort/s_proba


#val = calcul_puissance("stats_vent/stats_10kt_ter.csv", "polaires/BTS_11.csv", 11, True, "BTS", False)
technos = ["ADD", "BTS", "SS", "WISAMO", "ZEPHIRE"]
stats_transit = ["stats_vent/stats_10kt_transit_direct.csv","stats_vent/stats_transit_1WP.csv","stats_vent/stats_transit_2WP.csv"]
stats_fishing = ["stats_vent/stats_route_fishing_direct.csv","stats_vent/stats_route_fishing_1WP.csv","stats_vent/stats_route_fishing_2WP.csv"]
#stats_fishing = ["stats_vent/stats_route_fishing1.csv", "stats_vent/stats_route_fishing2.csv"]
#affichage_comparaison_route(stats)
#test = comparaison_route(11, technos, stats_fishing)
#comparaison_techno(11, technos, "stats_vent/stats_10kt_ter.csv", True)
#lecture_vent("stats_vent/stats_10kt_ter.csv",11,True,True)
#test = économie_puissance("puissances/sans_voile.csv", "puissances/avec_ADD.csv", "stats_vent/stats_10kt_bis.csv", 11,
 #                        20, True)
#test = lecture_vent("stats_vent/stats_10kt_ter.csv", 11)
# conversion_donnée(test, 11)
# print(test)
# lecture_vent("stats_vent/stats_10kt_bis.csv",True)
#test_bis = lecture_pol("polaires/ADD_11.csv")
#affichage_des_polaires(test_bis)
comparaison_route(9.5, technos, stats_transit)
# comparaison_route(11, technos, stats_fishing)
#affichage_comparaison_route(stats,11,False)
#affichage_comparaison_route(stats,11,True)
#test=comparaison_techno_aller_retour(11,technos,"stats_vent/stats_10kt.csv")
#print(test,"test")

#test=calc_puissance_pour_vitesse("stats_vent/stats_10kt_transit_direct.csv","polaires/ADD_11.csv",10,11)
# print(test)

lecture_vent("stats_vent/stats_10kt_transit_direct.csv",11,True,True)
# affichage_comparaison_route(stats_fishing,11,False)
# affichage_comparaison_route(stats_fishing,11,True)
