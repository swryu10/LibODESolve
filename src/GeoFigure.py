import os
import plotly.graph_objects as go
import GeoIATA

dir_build = '.'
name_exec = 'proj_GeodesicSphere.exec'

list_iata = []
list_city = []
list_lat = []
list_lon = []
fig = go.Figure(go.Scattergeo())

def reset() :
    list_iata = []
    list_city = []
    list_lat = []
    list_lon = []

    return
# end of function reset

def set_origin(code_iata_ini) :
    geo_loc_ini = GeoIATA.GeoLocation()
    geo_loc_ini.import_airport(code_iata_ini)
    if not geo_loc_ini.get_found() :
        return

    list_iata.append(code_iata_ini)
    str_iata_city = code_iata_ini + ' - ' + \
                    geo_loc_ini.get_city() + ', ' + \
                    geo_loc_ini.get_country()
    list_city.append(str_iata_city)
    list_lat.append(geo_loc_ini.get_latitude())
    list_lon.append(geo_loc_ini.get_longitude())

    return
# end of function set_origin

def add_destination(code_iata_add) :
    i_org = len(list_iata) - 1
    if i_org < 0 :
        return

    geo_loc_org = GeoIATA.GeoLocation()
    geo_loc_org.import_airport(list_iata[i_org])
    if not geo_loc_org.get_found() :
        return

    geo_loc_dst = GeoIATA.GeoLocation()
    geo_loc_dst.import_airport(code_iata_add)
    if not geo_loc_dst.get_found() :
        return

    i_dst = i_org + 1
    list_iata.append(code_iata_add)
    str_iata_city = code_iata_add + ' - ' + \
                    geo_loc_dst.get_city() + ', ' + \
                    geo_loc_dst.get_country()
    list_city.append(str_iata_city)
    list_lat.append(geo_loc_dst.get_latitude())
    list_lon.append(geo_loc_dst.get_longitude())

    str_command = dir_build + '/' + \
                  name_exec + ' ' + \
                  list_iata[i_org] + ' ' + \
                  list_iata[i_dst]
    print(str_command)
    val_return = os.system(str_command)

    name_file = 'geodesic_' + \
                list_iata[i_org] + '_' + \
                list_iata[i_dst] + '.txt'
    f_in = open(name_file, 'r')

    il = 0
    list_path_lat = []
    list_path_lon = []
    for line in f_in :
        list_word = line.split()
        if list_word[0][0] == '#' :
            continue

        if (il % 10 == 0) :
            list_path_lat.append(float(list_word[1]))
            list_path_lon.append(float(list_word[2]))

        il = il + 1

    fig.add_trace(
        go.Scattergeo(
            lat = list_path_lat,
            lon = list_path_lon,
            mode = 'lines',
            line = dict(width = 3, color = 'red')
        )
    )

    return
# end of function add_destination

def present() :
    fig.add_trace(
        go.Scattergeo(
            lat = list_lat,
            lon = list_lon,
            hoverinfo = 'text',
            text = list_city,
            mode = 'markers',
            marker = dict(
                size = 10,
                color = 'red',
                line = dict(
                    width = 3,
                    color = 'rgba(68, 68, 68, 0)'
                )
            )
        )
    )

    fig.update_geos(
        showcountries = True, countrycolor = "grey",
        showcoastlines = True, coastlinecolor = "black",
        showland = True, landcolor = "lightyellow",
        showocean = True, oceancolor = "lightskyblue",
        showlakes = True, lakecolor= "lightblue"
    )

    fig.update_geos(projection_type = "orthographic")
    fig.show()

    return
# end of function present
