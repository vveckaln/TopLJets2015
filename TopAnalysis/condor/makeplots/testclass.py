class Dog(object):

    def __init__(self, name):
        print "initiating class"
        self.name = name
        self.tricks = []    # creates a new empty list for each dog

    def add_trick(self, trick, var1, var2):
        print "var1 ", var1, "var2 ", var2
        self.tricks.append(trick)
        print self.tricks

d = Dog('Fido')
e = Dog('Buddy')
#d.add_trick()
d.add_trick('roll over', var2 = "a", var1 = "b")
try:
    #pass
    e.add_trickt('play dead')
except:
    print "some error occured"
print d.tricks
print e.tricks
