#!/usr/bin/python3
import socket
import datetime
import pytz
import numpy as np
import time
from astromath import AstroMath

from datetime import datetime, date, tzinfo
from math import pi


def sendPosition(ra_int, dec_int, status):
    size = 24
    rtype = 0
    datetime_now = int(time.time())
    buffer = bytearray()
    buffer.extend(size.to_bytes(2, byteorder="little"))
    buffer.extend(rtype.to_bytes(2, byteorder="little"))
    buffer.extend(datetime_now.to_bytes(8, byteorder="little"))
    buffer.extend(ra_int.to_bytes(4, byteorder="little"))
    buffer.extend(dec_int.to_bytes(4, byteorder="little"))
    buffer.extend(status.to_bytes(4, byteorder="little"))
    return bytes(buffer)


addr = socket.getaddrinfo("0.0.0.0", 10007)[0][-1]
# try:
#     rtc = RTC()
#     ntptime.settime()
# except:
#     print("Exception rtp")
s = socket.socket()
s.bind(addr)
s.listen(1)
print("listening on", addr)
while True:
    cl, addr = s.accept()
    print("client connected from", addr)
    while True:
        data = cl.recv(160)
        if data:
            size = int.from_bytes(data[0:2], byteorder="little")
            rtype = int.from_bytes(data[2:4], byteorder="little")
            if rtype == 0:
                client_micros = int.from_bytes(data[4:12], byteorder="little")
                client_time_s = client_micros / 1000000.0
                ra_int = int.from_bytes(data[12:16], byteorder="little")
                dec_int = int.from_bytes(data[16:20], byteorder="little", signed=True)
                ra_s = ra_int * (12 * 60 * 60 / 0x80000000)
                ra_rad = ra_int * (pi / 0x80000000)
                dec_deg = dec_int * (180 / 0x80000000)
                dec_rad = dec_int * (pi / 0x80000000)

                # print("length", size)
                # print("rtype", rtype)
                # print("client_micros", client_micros)
                date_tm = datetime.fromtimestamp(client_time_s).strftime("%H:%M:%S")
                # print("Client date time %s" % (date_tm))
                # print("Declination ", dec_deg)

                # print(
                #     "Rektaszension ",
                #     datetime.fromtimestamp(ra_s, tz=pytz.utc).strftime("%H:%M:%S.%f"),
                # )
                # EE =  –0,2317
                # x = print_ra + timedelta(seconds=3)
                # GAST0 = print_ra +EE
                # alpha = ra_int * (180 / 2147483648)
                longitude = 9.4542218
                latitude = 47.4926525
                am = AstroMath(longitude, latitude)
                now = datetime.utcnow()
                year = now.year
                month = now.month
                day = now.day
                utc = now.hour + now.minute / 60 + now.second / 3600
                year_dec = year + month / 12 + day / 365
                trans = am.Supplement(year_dec, 2000.0)
                ra_rad_trans, dec_rad_trans = am.Transform(ra_rad, dec_rad, trans)
                ra_s_trans = ra_rad_trans * 12 * 60 * 60 / pi
                dec_deg_trans = dec_rad_trans * 180 / pi

                latitude_rad = latitude * pi / 180
                height = 400
                gmst = am.SiderialTime(year, month, day, utc, 0)
                lmst = am.SiderialTime(year, month, day, utc, am.longitude)
                gmst_s = gmst * 60 * 60
                lmst_s = lmst * 60 * 60
                # print("Greenwich Mean Siderial Time (GMST): ", gmst)
                # print("Local Mean Siderial Time (LMST)    : ", lmst)
                # print(
                #     "Greenwich Mean Siderial Time (GMST): ",
                #     datetime.fromtimestamp(gmst_s, tz=pytz.utc).strftime("%H:%M:%S.%f"),
                # )

                print(
                    "Local Mean Siderial Time (LMST)    : ",
                    datetime.fromtimestamp(lmst_s, tz=pytz.utc).strftime("%H:%M:%S.%f"),
                )
                # Ruhende äquatoriale (δ,τ) in rotierende äquatoriale Koordinaten (δ,α) und umgekehrt
                # theta 	= Sternzeit am Ort der Beobachtung
                # tau 	    = Stundenwinkel
                # alpha 	= Rektaszension
                # delta 	= Deklination
                # Die Deklination δ bleibt unverändert.

                # alpha =theta -tau
                alpha = lmst_s - ra_s
                alpha_trans = lmst_s - ra_s_trans
                alpha_trans_rad = (alpha_trans * pi) / (12 * 60 * 60)
                am.PrecessToEOD(2000, ra_rad, dec_rad)
                print(
                    "RA",
                    datetime.fromtimestamp(ra_s, tz=pytz.utc).strftime("%H:%M:%S.%f"),
                    "Stundenwinkel",
                    datetime.fromtimestamp(alpha, tz=pytz.utc).strftime("%H:%M:%S.%f"),
                )
                print(
                    "RA Date",
                    datetime.fromtimestamp(ra_s_trans, tz=pytz.utc).strftime(
                        "%H:%M:%S.%f"
                    ),
                    "Stundenwinkel Date",
                    datetime.fromtimestamp(alpha_trans, tz=pytz.utc).strftime(
                        "%H:%M:%S.%f"
                    ),
                )
                print("Declination ", am.gms(dec_deg))
                print("Declination Date", am.gms(dec_deg_trans))
                az_rad, alt_rad = am.HA_DEC_to_AZ_ALT(
                    alpha_trans_rad, dec_rad_trans, latitude_rad
                )
                az = az_rad * 180 / pi
                alt = alt_rad * 180 / pi
                print("Azimut ", am.gms(az))
                print("Altitude", am.gms(alt))
                status = 0
                buff = sendPosition(ra_int, dec_int, status)
                # print("Buffer", buff)
                cl.send(buff)
                # tau =theta -alpha

            else:
                print("ops")
                break

    # cl.send(response)
    cl.close()
