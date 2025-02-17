import airportsdata
dic_airports = airportsdata.load("IATA")

class GeoLocation :
    def __init__(self) :
        self.code_iata = "NONE"
        self.found = False
        self.city = "NONE"
        self.country = "NONE"
        self.latitude = 0.0
        self.longitude = 0.0
    def import_airport(self, code_in_iata) :
        self.code_iata = code_in_iata
        self.found = self.code_iata in dic_airports
        if self.found :
            self.city = dic_airports[code_in_iata]["city"]
            self.country = dic_airports[code_in_iata]["country"]
            self.latitude = dic_airports[code_in_iata]["lat"]
            self.longitude = dic_airports[code_in_iata]["lon"]
    def get_found(self) :
        return self.found
    def get_city(self) :
        return self.city
    def get_country(self) :
        return self.country
    def get_latitude(self) :
        return self.latitude
    def get_longitude(self) :
        return self.longitude
