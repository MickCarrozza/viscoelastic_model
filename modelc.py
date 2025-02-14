class ModelC:
#   constructor
    def __init__(self,nr,lamb,Gmod,*nonlin):
        self.nr = nr
        self.lamb = lamb
        self.Gmod = Gmod

        match nr:
            case 1:
                self.__name = "Maxwell"
            case 2:
                self.__name = "Giesekus"
                self.alpha = nonlin[0]
                
        self.__created = True

#   function to check if created
    @property
    def is_created(self):
        return self.__created

#   function for model name
    @property
    def NAME(self):
        return self.__name

#   function for use with print(self)
    def __str__(self):
        return (f"{self.__name} model created with parameters:\n"
          f"lambda = {self.lamb}\n"
          f"G = {self.Gmod}\n")
