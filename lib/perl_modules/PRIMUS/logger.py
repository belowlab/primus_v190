'''-----------------------------------------
Adapted from versions by @eggins & @zruss11 
on github, respectively
--------------------------------------------
How to import this script into another file:
 	> from logger import logger, notify
 	> log = logger().log
* need to install dependencies first *
-----------------------------------------'''

import time, sys, os, platform
from datetime import datetime

class logger:

	def __init__(self):
		self.colors = {
			"error" 		: "\033[91m",
			"success" 		: "\033[92m",
			"info" 			: "\033[96m",
			"debug" 		: "\033[95m",
			"yellow" 		: "\033[93m",
			"lightpurple" 	: "\033[94m",
			"purp"			: "\033[94m",
			"lightgray" 	: "\033[97m",
			"clear"			: "\033[00m"
		}

	def log(self, message = "", color = "", script = "", file = "", procID = "", shown = True, showtime = True, nocolor = ""):
		
		currentTime = time.strftime("%H:%M:%S")

		if script and procID:	# for writing output to file in same format as API code
			timestring2 = '[%s] [%s] [%s] ' % (str(datetime.now()), "{: ^12}".format(script), procID)	
		else:
			timestring2 = '[%s] [%s] ' % (str(datetime.now()), "{: ^12}".format(script))

		try:
			colorString = self.colors[color]
		except:
			colorString = ""

		if showtime and script:
			timestring = "[%s] [%s] " % (currentTime, "{: ^12}".format(script))
		elif showtime:
			timestring = "[%s] " % currentTime			
		else:
			timestring = ""

		messageString = str(message) + self.colors['clear']
		noColorString = str(message)
		if nocolor:
			messageString += "%s" % nocolor 
			noColorString += "%s" % nocolor 

		finalString = "%s%s%s\n" % (timestring, colorString, str(messageString))
		noColorFinalString = "%s%s\n" % (timestring2, str(noColorString))

		sys.stdout.write(finalString)
		sys.stdout.flush()

		if file:
			with open(file, "a") as f:
				f.write(noColorFinalString)

# ------------------------------------------------------------------------------------
