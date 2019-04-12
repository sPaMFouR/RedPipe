
# # ------------------------------------------------------------------------------------------------------------------- #
# # Plot Light Curve
# # ------------------------------------------------------------------------------------------------------------------- #
# list_offset_mag_u = [float(value) + 1.5 for value in list_sorted_mag_u]
# list_offset_mag_b = [float(value) + 0.0 for value in list_sorted_mag_b]
# list_offset_mag_v = [float(value) - 0.8 for value in list_sorted_mag_v]
# list_offset_mag_r = [float(value) - 1.6 for value in list_sorted_mag_r]
# list_offset_mag_i = [float(value) - 2.5 for value in list_sorted_mag_i]
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Fit Spline To The Light Curves
# # ------------------------------------------------------------------------------------------------------------------- #
# plt.plot(list_offset_JD_u, list_sorted_mag_u, 'ko', markersize=8, label='$U$')
# plt.errorbar(list_offset_JD_u, list_sorted_mag_u, yerr=list_sorted_err_u, capsize=4, ls='none', color='k')
#
# plt.plot(list_offset_JD_b, list_sorted_mag_b, 'bD', markersize=8, label='$B$')
# plt.errorbar(list_offset_JD_b, list_sorted_mag_b, yerr=list_sorted_err_b, capsize=4, ls='none', color='k')
#
# plt.plot(list_offset_JD_v, list_sorted_mag_v, 'gs', markersize=8, label='$V$')
# plt.errorbar(list_offset_JD_v, list_sorted_mag_v, yerr=list_sorted_err_v, capsize=4, ls='none', color='k')
#
# plt.plot(list_offset_JD_r, list_sorted_mag_r, 'r^', markersize=8, label='$R$')
# plt.errorbar(list_offset_JD_r, list_sorted_mag_r, yerr=list_sorted_err_r, capsize=4, ls='none', color='k')
#
# plt.plot(list_offset_JD_i, list_sorted_mag_i, 'mp', markersize=8, label='$I$')
# plt.errorbar(list_offset_JD_i, list_sorted_mag_i, yerr=list_sorted_err_i, capsize=4, ls='none', color='k')
#
#
# def fit_spline1d(list_days, list_mag):
#     xaxis = np.linspace(min(list_days), max(list_days), 1000)
#     f1 = interpolate.interp1d(list_days, list_mag, kind='cubic')
#     plt.plot(xaxis, f1(xaxis), 'k-')
#
#
# def fit_curve_3d(list_days, list_mag):
#     xnew = np.linspace(min(list_days) - 2, max(list_days) + 2, 1000)
#     spl = interpolate.UnivariateSpline(list_days, list_mag, k=3, s=0.04)
#     plt.plot(xnew, spl(xnew), 'r-')
#
#
# def fit_curve_2d(list_days, list_mag):
#     xaxis = np.linspace(min(list_days), max(list_days), 1000)
#     spline = interpolate.UnivariateSpline(list_days, list_mag, k=2, s=0.1)
#     plt.plot(xaxis, spline(xaxis), 'k-')
#
#
# def fit_curve_1d(list_days, list_mag):
#     xaxis = np.linspace(min(list_days), max(list_days), 1000)
#     spline = interpolate.UnivariateSpline(list_days, list_mag, k=1, s=0.01)
#     plt.plot(xaxis, spline(xaxis), 'g-')
#
#
# def delete_element(list_of_lists, list_elements):
#     for element in list_elements:
#         for i in range(0, len(list_of_lists)):
#             list_of_lists[i] = [x for x in list_of_lists[i] if list_of_lists[i].index(x) != int(element)]
#     return list_of_lists
#
#
# [list_offset_JD_u, list_sorted_mag_u, list_sorted_err_u] = delete_element(
#     [list_offset_JD_u, list_sorted_mag_u, list_sorted_err_u], list_elements=[len(list_offset_JD_u) - 2])
#
# [list_offset_JD_b, list_sorted_mag_b, list_sorted_err_b] = delete_element(
#     [list_offset_JD_b, list_sorted_mag_b, list_sorted_err_b], list_elements=[len(list_offset_JD_b) - 1])
#
# [list_offset_JD_v, list_sorted_mag_v, list_sorted_err_v] = delete_element(
#     [list_offset_JD_v, list_sorted_mag_v, list_sorted_err_v], list_elements=[len(list_offset_JD_r) - 4])
#
# [list_offset_JD_r, list_sorted_mag_r, list_sorted_err_r] = delete_element(
#     [list_offset_JD_r, list_sorted_mag_r, list_sorted_err_r], list_elements=[len(list_offset_JD_r) - 2, 10])
#
# [list_offset_JD_i, list_sorted_mag_i, list_sorted_err_i] = delete_element(
#     [list_offset_JD_i, list_sorted_mag_i, list_sorted_err_i], list_elements=[20])
#
# # list_list_JD = [list_sorted_JD_u, list_sorted_JD_b, list_sorted_JD_v, list_sorted_JD_r, list_sorted_JD_i]
# list_list_JD = [list_offset_JD_u, list_offset_JD_b, list_offset_JD_v, list_offset_JD_r, list_offset_JD_i]
# list_list_mag = [list_sorted_mag_u, list_sorted_mag_b, list_sorted_mag_v, list_sorted_mag_r, list_sorted_mag_i]
#
# for index in range(0, len(list_list_JD)):
#     # fit_spline1d(list_list_JD[index], list_list_mag[index])
#     fit_curve_3d(list_list_JD[index], list_list_mag[index])
#     fit_curve_2d(list_list_JD[index], list_list_mag[index])
#     fit_curve_1d(list_list_JD[index], list_list_mag[index])
#
# plt.xlabel('$\Delta$t (Days Since $V$-Band Maximum)', fontsize=18)
# plt.ylabel('Apparent Magnitudes', fontsize=18)
# # plt.xlim(-10, 200)
# plt.ylim(22.0, 14.0)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# plt.legend(numpoints=1, fontsize=15)
# plt.title("Broadband Photometric Light Curves", fontsize=20)
# plt.grid()
# plt.savefig('fit_light_curve.png')
# plt.show()
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Print Interpolated Supernova Magnitudes In OUTPUT Files For Different Bands
# # ------------------------------------------------------------------------------------------------------------------- #
# def write_interpolated_sne_mag(list_jd, list_mag, output_file_name):
#     with open(str(output_file_name), "w") as f:
#         f.write(str('Julian_Day') + " " * 6 + str('JD_offset') + " " * 2 + str('SNe_Mag') + "\n")
#         jd_array = np.arange(int(min(list_jd)), int(max(list_jd)), 0.5)
#         spl = interpolate.UnivariateSpline(list_jd, list_mag, k=2, s=0.1)
#         for i in range(0, len(jd_array)):
#             f.write(str("%9.3f" % (float(jd_array[i]) + float(JD_offset))) + "     " +
#                     str("%7.3f" % (float(jd_array[i]))) + " " * 5 +
#                     str("%6.3f" % float(spl(jd_array[i]))) + "\n")
#
#
# write_interpolated_sne_mag(list_offset_JD_u, list_sorted_mag_u, output_file_name='OUTPUT_interp_snmag_u')
# write_interpolated_sne_mag(list_offset_JD_b, list_sorted_mag_b, output_file_name='OUTPUT_interp_snmag_b')
# write_interpolated_sne_mag(list_offset_JD_v, list_sorted_mag_v, output_file_name='OUTPUT_interp_snmag_v')
# write_interpolated_sne_mag(list_offset_JD_r, list_sorted_mag_r, output_file_name='OUTPUT_interp_snmag_r')
# write_interpolated_sne_mag(list_offset_JD_i, list_sorted_mag_i, output_file_name='OUTPUT_interp_snmag_i')
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Print Interpolated Supernova Magnitudes In A Single File
# # ------------------------------------------------------------------------------------------------------------------- #
# # with open('OUTPUT_interp_snmag', "w") as f:
# #     f.write(str('Julian_Day') + " " * 6 + str('JD_offset') + " " * 2 + str('SNe_Mag') + "\n")
# #     list_spl = []
# #     for list_jd in list_list_JD:
# #         for list_mag in list_list_mag:
# #             list_spl.append(interpolate.UnivariateSpline(list_jd, list_mag, k=3, s=0.04))
# #
# #             jd_array = np.arange(int(min(list_jd)), int(max(list_jd)), 0.5)
# #     for i in range(0, len(jd_array)):
# #         for list_mag
# #         f.write(str("%9.3f" % float(jd_array[i])) + "     " +
# #                 str("%7.3f" % (float(jd_array[i]) - float(JD_offset))) + " " * 5 +
# #                 str("%6.3f" % float(spl(jd_array[i]))) + "\n")
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # # ------------------------------------------------------------------------------------------------------------------- #
# # # Calculate Flux Values In Individual Bands Day-Wise
# # # ------------------------------------------------------------------------------------------------------------------- #
# # # mag_lambda = -2.5log(flux_lambda) - 21.100 - zp(f_lambda)
# # # Zero Point        U       B       V       R       I
# # # zp(f_lambda)   -0.152  -0.602   0.000   0.555   1.271
# # zp_u = -0.152
# # zp_b = -0.602
# # zp_v = 0.000
# # zp_r = 0.555
# # zp_i = 1.271
# # A_lambda_u = 0.303
# # A_lambda_b = 0.254
# # A_lambda_v = 0.192
# # A_lambda_r = 0.152
# # A_lambda_i = 0.105
# # list_zp = [zp_u, zp_b, zp_v, zp_r, zp_i]
# #
# #
# # def calculate_flux(mag_supernova_day):
# #     flux_u = 'INDEF'
# #     flux_b = 'INDEF'
# #     flux_v = 'INDEF'
# #     flux_r = 'INDEF'
# #     flux_i = 'INDEF'
# #
# #     if mag_supernova_day[0] != 'INDEF':
# #         flux_u = 10 ** (-(float(mag_supernova_day[0]) - A_lambda_u + zp_u + 21.100) / 2.5)
# #     if mag_supernova_day[1] != 'INDEF':
# #         flux_b = 10 ** (-(float(mag_supernova_day[1]) - A_lambda_b + zp_b + 21.100) / 2.5)
# #     if mag_supernova_day[2] != 'INDEF':
# #         flux_v = 10 ** (-(float(mag_supernova_day[2]) - A_lambda_v + zp_v + 21.100) / 2.5)
# #     if mag_supernova_day[3] != 'INDEF':
# #         flux_r = 10 ** (-(float(mag_supernova_day[3]) - A_lambda_r + zp_r + 21.100) / 2.5)
# #     if mag_supernova_day[4] != 'INDEF':
# #         flux_i = 10 ** (-(float(mag_supernova_day[4]) - A_lambda_i + zp_i + 21.100) / 2.5)
# #     return [flux_u, flux_b, flux_v, flux_r, flux_i]
# #
# # flux_supernova = []
# # for index in range(0, len(mag_supernova)):
# #     flux_supernova.append(calculate_flux(mag_supernova[index]))
# #
# #
# # with open('output_flux_values' , 'w') as f:
# #     for i in range(0, len(flux_supernova)):
# #         for j in range(0, len(flux_supernova[i])):
# #             f.write(list_JD[j] + ' ' + str(flux_supernova[i][j]) + '\n')
# #         f.write('\n')
# #
# # # left_limit = 300
# # # right_limit = 1100
# # # list_band_centers = [370, 420, 530, 600, 800]
# # # list_band_centers.sort()
# # # test_y = [flux_supernova[0][0], flux_supernova[0][1], flux_supernova[0][2], flux_supernova[0][4]]
# # # test_x = [list_band_centers[0], list_band_centers[1], list_band_centers[2], list_band_centers[4]]
# # #
# # # for i in range(0, len(flux_supernova)):
# # #     print flux_supernova[i]
# # # plt.plot(list_band_centers, flux_supernova[2],'o')
# # # spline = interpolate.UnivariateSpline(list_band_centers, flux_supernova[0], bbox=[300, 800], k=4)
# # # spline2 = interpolate.interp1d(list_band_centers, flux_supernova[0], kind='linear')
# # # spline3 = interpolate.splrep(test_x, test_y, s=2, k=3)
# # # xaxis = np.linspace(200, 1000, 1000)
# # # ynew = interpolate.splev(xaxis, spline3, der=0)
# # # plt.plot(xaxis, ynew)
# # # plt.ylim(-0.1e-15, 2e-15)
# # # plt.grid()
# # # plt.show()
# # # ------------------------------------------------------------------------------------------------------------------- #
# #
# # for i in range(0, len(list_dates)):
# #     with open("output_flux_" + str(list_dates[i]), 'w') as f:
# #         f.write(str(310) + " " + str(920) + " " + str(100) + "\n")
# #         if flux_supernova[i][0] != 'INDEF':
# #             f.write(str(list_band_centers[0]) + " " + str(flux_supernova[i][0]) + "\n")
# #         if flux_supernova[i][1] != 'INDEF':
# #             f.write(str(list_band_centers[1]) + " " + str(flux_supernova[i][1]) + "\n")
# #         if flux_supernova[i][2] != 'INDEF':
# #             f.write(str(list_band_centers[2]) + " " + str(flux_supernova[i][2]) + "\n")
# #         if flux_supernova[i][3] != 'INDEF':
# #             f.write(str(list_band_centers[3]) + " " + str(flux_supernova[i][3]) + "\n")
# #         if flux_supernova[i][4] != 'INDEF':
# #             f.write(str(list_band_centers[4]) + " " + str(flux_supernova[i][4]) + "\n")
# #
# # # print list_fits_v[14:17]
# # # print list_fits_b[6:10]
# # #
# # # print mag_supernova_v[14:17]
# # # print mag_supernova_b[6:10]
# # # trapz(a[:, 1], x=a[:, 0])
# # bol_flux_dec19 = 2.672e-13
# # bol_flux_dec26 = 2.287e-13
# # bol_flux_jan03 = 1.973e-13
# # bol_flux_jan25 = 1.794e-13
# # bol_flux_jul25 = 5.446e-12
# # bol_flux_jul16 = 7.703e-12
# # bol_flux_175 = 2.25e-13
# # dist = 45 * 3.086e+24
# # z = 0.010424
# #
# # test_lum1 = 4 * np.pi * dist**2 * bol_flux_dec19
# # test_lum2 = 4 * np.pi * dist**2 * bol_flux_jan25
# # test_mag1 = -2.5 * math.log10(test_lum1 / 3.839e+33) + 4.72
# # test_mag2 = -2.5 * math.log10(test_lum2 / 3.839e+33) + 4.72
# # bol_lum_175 = bol_flux_175 * 4 * np.pi * (dist ** 2)
# # print "bol_lum_dec19 = ", test_lum1
# # print "bol_lum_jan25 = ", test_lum2
# # print "\n"
# # print -0.4 * math.log10(test_lum1)
# #
# # print "bol_mag_dec19 = ", test_mag1
# # print "bol_mag_jan25 = ", test_mag2
# # print "\n"
# #
# # m_ni2 = (bol_lum_175 * 0.07 / 9.92e41) / (math.exp(-175/111.4) - math.exp(-175/8.8))
# # nickel_mass = 7.866e-44 * bol_lum_175 * math.exp(((175 * (1 + z) - 6.1) / 111.26))
# #
# # print "jerkstrand_Ni = ", m_ni2
# # print "hamuy_ni =", nickel_mass
# #
# # y = [bol_flux_dec19, bol_flux_dec26, bol_flux_jan03, bol_flux_jan25]
# # x = [list_JD_v[10] + 11, list_JD_v[11] + 11, list_JD_v[12] + 11, list_JD_v[13] + 11]
# # # plt.plot(x, y, 'o')
# # # y2 = interpolate.interp1d(x, y, kind='linear')
# # # plt.plot(x, y2(x), '-')
# # # # plt.show()
# #
# # peak_bolometric_lum = 4 * np.pi * dist**2 * bol_flux_jul16
# # peak_bolometric_mag = -2.5 * math.log10(peak_bolometric_lum / 3.839e+33) + 4.72
# # print -0.4 * math.log10(peak_bolometric_lum)
# # print "peak bolometric lum = ", peak_bolometric_lum
# # print "peak bolometric mag =", peak_bolometric_mag
# #
# # tail_lum = bol_lum_175
# # print "tail_lum=", tail_lum
# # tail_mag_et = 1e40 * 10**(0.8)
# # print tail_lum/tail_mag_et * 0.048
# #
# # lum_175 = 10**(-((20.0 - 0.26 - 33.2) - 4.72)/2.5) * 3.839e33
# # print "tail_lum=", lum_175
# #
# # # v = 2000 kms-1
# # # velociy around v-50
# # # spiro et al 2014 mnras
# # # Intro - Finding chart
# # #
# # # Light Curve
# # # ColorCurves
# # # fitting the bolometric flux
# # # Estimation of ni mass
# # # Nickel Mass - duration of plateau light curve
# # # v -band steepness - slope of 4 points
# # # nickel mass - radioactive tail
# # # observations - type IIP supernova
# # # duartion of the plateau
# # #
# # # vinko et al
# # #
# # # absolute v magntitude = 17 - 0.2 - 33.2
# # # low luminosity type-IIP
# # # 1.25 V-I
# # # 0.60 V-R
# # # 1.00 B-V
# #
# # list_steep = list_JD_v[23:27]
# # mag_steep = mag_supernova_v[23:27]
# #
# # plt.plot(list_steep, mag_steep)
# #
# # slope = (mag_steep[3] - mag_steep[0]) / (list_steep[3] - list_steep[0])
# # # print list_steep
# # # print mag_steep
# # print slope
# # M_ni_2 = 10 ** (-6.2295 * slope - 0.8417)
# # print M_ni_2
#
# # s1rev = interpolate.InterpolatedUnivariateSpline(list_days, list_mag)
# # s3 = interpolate.interp1d(list_days, list_mag, kind='cubic')
# # rbf = interpolate.Rbf(list_days, list_mag)
# # tck = interpolate.splrep(list_days, list_mag, s=1)
# # ynew = interpolate.splev(xaxis, tck, der=0)
# # x = np.arange(0, 2*np.pi+np.pi/4, 2*np.pi/8)
# # y = np.sin(x)
# # s = interpolate.InterpolatedUnivariateSpline(x, y)
# # xnew = np.arange(0, 2*np.pi, np.pi/50)
# # ynew = s(xnew)
# #
# # plt.figure()
# # plt.plot(x, y, 'x', xnew, ynew, xnew, np.sin(xnew), x, y, 'b')
# # plt.legend(['Linear', 'InterpolatedUnivariateSpline', 'True'])
# # plt.axis([-0.05, 6.33, -1.05, 1.05])
# # plt.title('InterpolatedUnivariateSpline')
# # plt.show()
#
# # m_star - m_sun = -2.5*math.log10(f_star/f_sun)
# # f_star = f_sun * 10 ** (-0.4 * (m_star - m_sun))
# # f_star = 3.839e33 * 10 ** (-0.4 * (m - 4.72))
