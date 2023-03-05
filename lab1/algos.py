import math

def euler(x, y, h, n, f):
    ans = []
    for i in range(n):
        try:
            y += h * f(x, y)
            ans.append(y)
            x += h
        except OverflowError:
            ans.append("Over")
    return ans

#---------------------------------------------------#

def picar_apox_11(u):
    return 1 + u + u**3 / 3

def picar_apox_12(u):
    return picar_apox_11(u) + u**2 / 2 + u**4 / 12

def picar_apox_13(u):
    return picar_apox_12(u) + u**3 / 6 + u**5 / 60

def picar_apox_14(u):
    return picar_apox_13(u) + u**4 / 24 + u**6 / 360

#---------------------------------------------------#

def picar_apox_21(u):
    return 0.5 + u**2 / 2 + u**4 / 4

def picar_apox_22(u):
    return picar_apox_21(u) + u**4 / 4 + u**6 / 12

def picar_apox_23(u):
    return picar_apox_22(u) + u**6 / 12 + u**8 / 48

def picar_apox_24(u):
    return picar_apox_23(u) + u**8 / 48 + u**10 / 240

#---------------------------------------------------#

def picar_apox_31(u):
    return u**3 / 3

def picar_apox_32(u):
    return picar_apox_31(u) + u**7 / 63

def picar_apox_33(u):
    return picar_apox_32(u) + 2 * u**11 / 2079 + u**15 / 59535

def picar_apox_34(u):
    return picar_apox_33(u) + (2 / 93555)*u**15 + (2 / 3393495)*u**19 + \
        (2 / 2488563)*u** 19 + (2 / 86266215)*u**23 + (1 / 99411543)*u**23 + \
        (2 / 3341878155)*u**27 + (1 / 109876902975)*u**31

#---------------------------------------------------#

def picar(x_max, h, apox_fun, x, u):
    ans = []
    while abs(x) < abs(x_max):
        ans.append(u)
        x += h
        u = apox_fun(x)
    return ans