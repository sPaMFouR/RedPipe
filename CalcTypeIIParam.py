#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxx----------------------CALCULATE TYPE II SN PARAMS---------------------xxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import numpy as np
import pandas as pd
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Change The Observed Parameters Here
# ------------------------------------------------------------------------------------------------------------------- #
Rv = 3.1
redshift = 0.008059
NaID_MW = 0.44
NaID_MWErr = 0.08
NaID_host = 0.89
NaID_hostErr = 0.13
VIcolor = 0.831
VIcolorErr = 0.107
VIphase = 84.0
phase = 200
steepness = 0.121
deltat = 88
deltatErr = 5
vph = 3913
vphErr = 320
mv50 = -17.16
mv50Err = 0.37


EBV_mag = 0.21
EBV_Err = 0.11

H0 = 73.52
Ab = 4.06 * EBV_mag
AbErr = 4.06 * EBV_Err

Av = 3.04 * EBV_mag
AvErr = 3.04 * EBV_Err

Ai = 1.71 * EBV_mag
AiErr = 1.71 * EBV_Err

VIplus50 = 0.52
VIplus50Err = 0.02
VIminus30 = 0.64
VIminus30Err = 0.02

vplus50 = 4272
vplus50Err = 53
vminus30 = 3022
vminus30Err = 42

Bplus50 = 17.19
Bplus50Err = 0.01
Bminus30 = 17.48
Bminus30Err = 0.02

Vplus50 = 16.26
Vplus50Err = 0.01
Vminus30 = 16.29
Vminus30Err = 0.01

Iplus50 = 15.47
Iplus50Err = 0.02
Iminus30 = 15.37
Iminus30Err = 0.02

dict_bol = {198.7: {'BolLum': 1.047e41, 'BolErr': 1.504e40}, 212.8: {'BolLum': 9.694e40, 'BolErr': 1.393e40},
            237.8: {'BolLum': 8.064e40, 'BolErr': 1.158e40}, 240.8: {'BolLum': 7.012e40, 'BolErr': 1.008e40}}
data_bol = pd.DataFrame(dict_bol).T
data_bol.index.name = 'Phase'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Calculating Parameters
# ------------------------------------------------------------------------------------------------------------------- #

def hamuyni(lum, phase):
    return 7.866e-44 * lum * np.exp((phase * (1 + redshift) - 6.1) / 111.26)


def get_extbarbon(ew, ew_Err=0.0):
    ebv = 0.25 * ew
    ebv_Err = 0.25 * ew_Err
    print (r"E(B-V) [Barbon et al. (1990)] = {0:0.2f}+/-{1:0.2f}".format(ebv, ebv_Err))


def get_extturatto(ew, ew_Err=0.0):
    ebv = -0.01 + 0.16 * ew
    ebv_Err = 0.16 * ew_Err
    print (r"E(B-V) [Turatto et al. (2003)] = {0:0.2f}+/-{1:0.2f}".format(ebv, ebv_Err))


def get_extpoznanski(ew, ew_Err=0.0):
    logebv = 1.17 * ew - 1.85
    logebv_Err = (1.17 * ew_Err ** 2 + 0.08 ** 2) ** 0.5
    ebv = 10 ** logebv
    ebv_Err = 10 ** (logebv + logebv_Err) - (10 ** logebv)
    print (r"E(B-V) [Poznanski et al. (2012)] = {0:0.2f}+/-{1:0.2f}".format(ebv, ebv_Err))


def get_exthostcolmethod(color, colorErr):
    Av = 2.518 * (color - 0.656)
    AvErr = 2.518 * (colorErr ** 2 + 0.053 ** 2 + 0.059 ** 2) ** 0.5
    print (r"E(B-V) [Olivares et al. (2010)] = {0:0.2f}+/-{1:0.2f}".format(Av / Rv, AvErr / Rv))


def get_hamuyni(data_bol):
    lum = data_bol['BolLum'].mean()
    phase = np.mean(data_bol.index.values)
    lumErr1 = data_bol['BolLum'].std()
    lumErr2 = data_bol['BolErr'].apply(lambda x: x / data_bol.shape[0]).mean()
    nimass = hamuyni(lum, phase)
    niErr = hamuyni(lum + (lumErr1 ** 2 + lumErr2 ** 2) ** 0.5, phase) - nimass
    print "Mean Phase = {0:0.1f}".format(phase)
    print (r"Mean Tail Luminosity {0:0.3e}+/-{1:0.3e}".format(lum, (lumErr1 ** 2 + lumErr2 ** 2) ** 0.5))
    print (r"$^{56}Ni$ Mass [Hamuy (2003)]" + " = {0:0.3f}+/-{1:0.3f}".format(nimass, niErr))


def get_jerkni(lum, phase):
    nimass = (lum * 0.07 / 9.92e41) / (np.exp(- phase / 111.26) - np.exp(- phase / 8.8))
    print (r"$^{56}Ni$ Mass [Hamuy (2003)]" + " = {0:0.3f}".format(nimass))


def get_steepnessni(steepness):
    # lognimass = -3.0805 * steepness - 1.0861
    # logniErr = ((0.1224 * steepness) ** 2 + 0.0196 ** 2) ** 0.5
    lognimass = -3.5024 * steepness - 1.0167
    logniErr = ((0.0960 * steepness) ** 2 + 0.003 ** 2) ** 0.5
    nimass = 10 ** lognimass
    niErr = 10 ** (lognimass + logniErr) - 10 ** lognimass
    print (r"$^{56}Ni$ Mass [Singh et al. (2018)]" + " = {0:0.3f}+/-{1:0.3f}".format(nimass, niErr))


def get_litvinovaparams(mv, deltat, vph):
    vph = vph / 1000
    logE = 0.135 * mv + 2.34 * np.log10(deltat) + 3.13 * np.log10(vph) - 3.205
    logM = 0.234 * mv + 2.91 * np.log10(deltat) + 1.96 * np.log10(vph) - 1.829
    logR = -0.572 * mv - 1.07 * np.log10(deltat) - 2.74 * np.log10(vph) - 3.350

    print (r"Explosion Energy = {0:0.3f}".format(10 ** logE))
    print (r"Ejecta Mass = {0:0.3f}".format(10 ** logM))
    print (r"Progenitor Radius = {0:0.3f}".format(10 ** logR))


def get_scmhamuy(band, Alambda, AlambdaErr, vphot, vErr, mag, magErr):
    if band == 'V':
        logD = 0.2 * (mag - Alambda + 6.25 * np.log10(vphot / 5000.) + 1.46)
        logDErr = 0.2 * (magErr ** 2 + AlambdaErr ** 2 + (1.35 * np.log10(vphot / 5000.)) ** 2 +
                         (6.25 * vErr / vphot) ** 2 + 0.15 ** 2) ** 0.5
    elif band == 'I':
        logD = 0.2 * (mag - Alambda + 5.45 * np.log10(vphot / 5000.) + 1.92)
        logDErr = 0.2 * (magErr ** 2 + AlambdaErr ** 2 + (0.91 * np.log10(vphot / 5000.)) ** 2 +
                         (5.45 * vErr / vphot) ** 2 + 0.11 ** 2) ** 0.5
    else:
        print ("ERROR: SCM Technique Not Valid For Band '{0}'".format(band))
        sys.exit(1)

    D = (10 ** logD) / H0
    DErr = (10 ** (logD + logDErr)) / H0 - D

    print (r"Distance [Hamuy et al. (2004)] = {0:0.2f}+/-{1:0.2f}".format(D, DErr))
    print (r"$\mu$ = " + "{0:0.2f}+/-{1:0.2f}".format(5 * np.log10(D * 1e6) - 5, 5 * np.log10((D + DErr) / D)))

    return D, DErr


def get_scmnugent(band, vphot, vErr, mag, magErr, color, colorErr):
    if band == 'I':
        logD = mag + 6.69 * np.log10(vphot / 5000.) + 1.36 * (color - 0.53) + 17.49
        logDErr = (magErr ** 2 + (0.50 * np.log10(vphot / 5000.)) ** 2 + (6.69 * vErr / vphot) ** 2 +
                   (1.36 * colorErr) ** 2 + 0.08 ** 2) ** 0.5
    else:
        print ("ERROR: SCM Technique Not Valid For Band '{0}'".format(band))
        sys.exit(1)

    D = 10 ** ((logD + 5) / 5) / 1e6
    DErr = (10 ** ((logD + logDErr + 5) / 5) - 10 ** ((logD + 5) / 5)) / 1e6

    print (r"Distance [Nugent et al. (2006)] = {0:0.2f}+/-{1:0.2f}".format(D, DErr))
    print (r"$\mu$ = " + "{0:0.2f}+/-{1:0.2f}".format(logD, logDErr))

    return D, DErr


def get_scmpoznanski(band, vphot, vErr, mag, magErr, color, colorErr):
    if band == 'I':
        logD = 0.2 * (mag + 4.4 * np.log10(vphot / 5000.) - 0.80 * (color - 0.53) + 1.76)
        logDErr = 0.2 * (magErr ** 2 + (0.60 * np.log10(vphot / 5000.)) ** 2 + (4.4 * vErr / vphot) ** 2 +
                         (0.80 * colorErr) ** 2 + (0.3 * color) ** 2 + 0.05 ** 2) ** 0.5
    else:
        print ("ERROR: SCM Technique Not Valid For Band '{0}'".format(band))
        sys.exit(1)

    D = (10 ** logD) / H0
    DErr = (10 ** (logD + logDErr)) / H0 - D

    print (r"Distance [Poznanski et al. (2009)] = {0:0.2f}+/-{1:0.2f}".format(D, DErr))
    print (r"$\mu$ = " + "{0:0.2f}+/-{1:0.2f}".format(5 * np.log10(D * 1e6) - 5, 5 * np.log10((D + DErr) / D)))

    return D, DErr


def get_scmolivares(band, vphot, vErr, mag, magErr, color, colorErr):
    if band == 'B':
        logD = 0.2 * (mag + 3.50 * np.log10(vphot / 5000.) + 1.99 - 2.67 * color)
        logDErr = 0.2 * (magErr ** 2 + (0.30 * np.log10(vphot / 5000.)) ** 2 + (3.50 * vErr / vphot) ** 2 +
                         0.11 ** 2 + (0.13 * color) ** 2 + (2.67 * colorErr) ** 2) ** 0.5
    elif band == 'V':
        logD = 0.2 * (mag + 3.08 * np.log10(vphot / 5000.) + 2.38 - 1.67 * color)
        logDErr = 0.2 * (magErr ** 2 + (0.25 * np.log10(vphot / 5000.)) ** 2 + (3.08 * vErr / vphot) ** 2 +
                         0.09 ** 2 + (0.10 * color) ** 2 + (1.67 * colorErr) ** 2) ** 0.5
    elif band == 'I':
        logD = 0.2 * (mag + 2.62 * np.log10(vphot / 5000.) + 2.23 - 0.60 * color)
        logDErr = 0.2 * (magErr ** 2 + (0.21 * np.log10(vphot / 5000.)) ** 2 + (2.62 * vErr / vphot) ** 2 +
                         0.07 ** 2 + (0.09 * color) ** 2 + (0.60 * colorErr) ** 2) ** 0.5
    else:
        print ("ERROR: SCM Technique Not Valid For Band '{0}'".format(band))
        sys.exit(1)

    D = (10 ** logD) / H0
    DErr = (10 ** (logD + logDErr)) / H0 - D

    print (r"Distance [Olivares et al. (2010)] = {0:0.2f}+/-{1:0.2f}".format(D, DErr))
    print (r"$\mu$ = " + "{0:0.2f}+/-{1:0.2f}".format(5 * np.log10(D * 1e6) - 5, 5 * np.log10((D + DErr) / D)))

    return D, DErr

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculating Parameters
# ------------------------------------------------------------------------------------------------------------------- #

get_extbarbon(NaID_MW, NaID_MWErr)
get_extbarbon(NaID_host, NaID_hostErr)

get_extturatto(NaID_MW, NaID_MWErr)
get_extturatto(NaID_host, NaID_hostErr)

get_extpoznanski(NaID_MW, NaID_MWErr)
get_extpoznanski(NaID_host, NaID_hostErr)

get_exthostcolmethod(VIcolor, VIcolorErr)

get_hamuyni(data_bol)
get_steepnessni(steepness)

get_litvinovaparams(mv50, deltat, vph)

dist = np.zeros(7)
dist_err = np.zeros(7)

dist[0], dist_err[0] = get_scmhamuy('V', Av, AvErr, vplus50, vplus50Err, Vplus50, Vplus50Err)
dist[1], dist_err[1] = get_scmhamuy('I', Ai, AiErr, vplus50, vplus50Err, Iplus50, Iplus50Err)

dist[2], dist_err[2] = get_scmnugent('I', vplus50, vplus50Err, Iplus50, Iplus50Err, VIplus50, VIplus50Err)
dist[3], dist_err[3] = get_scmpoznanski('I', vplus50, vplus50Err, Iplus50, Iplus50Err, VIplus50, VIplus50Err)

dist[4], dist_err[4] = get_scmolivares('B', vminus30, vminus30Err, Bminus30, Bminus30Err, VIminus30, VIminus30Err)
dist[5], dist_err[5] = get_scmolivares('V', vminus30, vminus30Err, Vminus30, Vminus30Err, VIminus30, VIminus30Err)
dist[6], dist_err[6] = get_scmolivares('I', vminus30, vminus30Err, Iminus30, Iminus30Err, VIminus30, VIminus30Err)

print np.round(np.mean(dist), 2)
print np.round((np.std(dist) ** 2 + np.sum(np.square(dist_err)) / 7 ** 2) ** 0.5, 2)
# ------------------------------------------------------------------------------------------------------------------- #
