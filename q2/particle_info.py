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
            self.leaves = True

        def g(self):
            return 7 / 8 * 4 * 3

        def T(self):
            return 30 * 1000

    class bottom:
        def __init__(self):
            self.leaves = True

        def g(self):
            return 7 / 8 * 4 * 3

        def T(self):
            return 1 * 1000  # 1 GeV

    class charm:
        def __init__(self):
            self.leaves = True

        def g(self):
            return 7 / 8 * 4 * 3

        def T(self):
            return 200

    class strange:
        def __init__(self):
            self.leaves = True

        def g(self):
            return 7 / 8 * 4 * 3

        def T(self):
            return 200

    class up:
        def __init__(self):
            self.leaves = True

        def g(self):
            return 7 / 8 * 4 * 3

        def T(self):
            return 200

    class down:
        def __init__(self):
            self.leaves = True

        def g(self):
            return 7 / 8 * 4 * 3

        def T(self):
            return 200


class gluon:
    def __init__(self):
        self.leaves = True

    def g(self):
        return 8 * 2

    def T(self):
        return 200


class boson:
    def __init__(self):
        self.W = self.W()
        self.Z = self.Z()
        self.h = self.h()

    class W:
        def __init__(self):
            self.leaves = True

        def g(self):
            return 3 * 2

        def T(self):
            return 15 * 1000

    class Z:
        def __init__(self):
            self.leaves = True

        def g(self):
            return 3

        def T(self):
            return 15 * 1000

    class h:
        def __init__(self):
            self.leaves = True

        def g(self):
            return 1

        def T(self):
            return 15 * 1000


class X:
    def __init__(self):
        self.leaves = True

    def g(self):
        return 4

    def T(self):
        return 100 * 1000


class pion:
    def __init__(self):
        self.pi_pm = self.pi_pm()
        self.pi_0 = self.pi_0()

    class pi_pm:
        def __init__(self):
            pass

        def g(self):
            return 2

        def T_join(self):
            return 200

        def T_leave(self):
            return 50

    class pi_0:
        def __init__(self):
            pass

        def g(self):
            return 2

        def T_join(self):
            return 200

        def T_leave(self):
            return 50


class lepton:
    def __init__(self):
        self.e = self.e()
        self.muon = self.muon()
        self.tau = self.tau()

    class e:
        def __init__(self):
            self.leaves = True

        def g(self):
            return 7 / 8 * 4

        def T(self):
            return 0.25

    class muon:
        def __init__(self):
            self.leaves = True

        def g(self):
            return 7 / 8 * 4

        def T(self):
            return 50

    class tau:
        def __init__(self):
            self.leaves = True

        def g(self):
            return 7 / 8 * 4

        def T(self):
            return 200


class photon:
    def __init__(self):
        pass

    def g(self):
        return 2

    def T(self):
        pass


class neutrino:
    def __init__(self):
        pass

    def g(self):
        return 7 / 8 * 6 * (4 / 11) ** (4 / 3)

    def g_s(self):
        return 7 / 8 * 6 * (4 / 11)

    def T(self):
        pass
