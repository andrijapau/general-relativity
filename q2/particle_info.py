class quark:
    def __init__(self):
        self.top = self.top()
        self.bottom = self.bottom()
        self.charm = self.charm()
        self.strange = self.strange()
        self.up = self.up()
        self.down = self.down()

    class top:
        def __init__(self):
            pass
        def g(self):
            return 7/8 * 4 * 3
        def critical_temp(self):
            return 30 * 1000


    class bottom:
        def __init__(self):
            pass
        def g(self):
            return 7/8 * 4 * 3
        def critical_temp(self):
            return 1 * 1000 # 1 GeV


    class charm:
        def __init__(self):
            pass
        def g(self):
            return 7/8 * 4 * 3
        def critical_temp(self):
            return 200


    class strange:
        def __init__(self):
            pass
        def g(self):
            return 7/8 * 4 * 3
        def critical_temp(self):
            return 200


    class up:
        def __init__(self):
            pass
        def g(self):
            return 7/8 * 4 * 3
        def critical_temp(self):
            return 200


    class down:
        def __init__(self):
            pass
        def g(self):
            return 7/8 * 4 * 3
        def critical_temp(self):
            return 200 

class gluon:
    def __init__(self):
        pass
    def g(self):
        return 8 * 2
    def critical_temp(self):
        return 200 

class W:
    def __init__(self):
        pass
    def g(self):
        return 3 * 2
    def critical_temp(self):
        return 15 * 1000

class Z:
    def __init__(self):
        pass
    def g(self):
        return 3
    def critical_temp(self):
        return 15 * 1000 

class h:
    def __init__(self):
        pass
    def g(self):
        return 1
    def critical_temp(self):
        return 15 * 1000 

class X:
    def __init__(self):
        pass
    def g(self):
        return 4
    def critical_temp(self):
        return 100 * 1000 

class pi_pm:
    def __init__(self):
        pass
    def g(self):
        return 2
    def critical_temp(self):
        return 200 

class pi_0:
    def __init__(self):
        pass
    def g(self):
        return 2
    def critical_temp(self):
        return 200 

class lepton:
    def __init__(self):
        self.e = self.e()
        self.muon = self.muon()
        self.tau = self.tau()

    class e:
        def __init__(self):
            pass
        def g(self):
            return 7/8 * 4
        def critical_temp(self):
            return 250 / 1000 

    class muon:
        def __init__(self):
            pass
        def g(self):
            return 7/8 * 4
        def critical_temp(self):
            return 50  

    class tau:
        def __init__(self):
            pass
        def g(self):
            return 7/8 * 4
        def critical_temp(self):
            pass

class photon:
    def __init__(self):
        pass
    def g(self):
        return 2
    def critical_temp(self):
        pass

class neutrino:
    def __init__(self):
        pass
    def g(self):
        return 7/8 * 2
    def critical_temp(self):
        pass

quark = quark()
gluon = gluon()
W_boson = W()
Z_boson = Z()
h_boson = h()
pi_pm = pi_pm()
pi_0 = pi_0()
lepton = lepton()
photon = photon()
neutrino = neutrino()
