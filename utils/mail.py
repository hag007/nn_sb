"""The first step is to create an SMTP object, each object is used for connection
with one server."""

import smtplib
server = smtplib.SMTP('smtp.gmail.com', 587)
server.ehlo()
server.starttls()
server.login("hagai.levi.007", "b6437355")
server.ehlo()

#Send the mail
msg = "Hello!" # The /n separates the message from the headers
server.sendmail("hagai.levi.007@gmail.com", "hagai.levi.007@gmail.com", msg)