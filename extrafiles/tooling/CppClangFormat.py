#!/usr/bin/env python3
"""
This tool iterates over all C++ files and applies clang-format
to them. Source files need to opt-in by containing the term
"networkit-format" (preferably in a comment near the header).

If used in read-only mode it returns a negative exit code if
a change is necessary.
"""
import filecmp
import nktooling as nkt
import tempfile
import sys
import shutil
import subprocess
import os

def getEnv():
	"""
	Returns the current environment variables with 'LANG' set to 'C' to ensure
	the output to be in English.
	"""
	env = dict(os.environ)
	env['LANG'] = 'C'
	return env

def isSupported(cmd):
	""" Checks if the given clang-format command is available and if its version is recent enough. """
	if shutil.which(cmd) is None:
		return False
	# Read major revision number from "clang-format version XX.XX.XX ... "
	version = str(subprocess.check_output([cmd, "--version"],
		universal_newlines=True, env=getEnv())).strip().split()[-1].split('.')[0]
	return int(version) >= nkt.MIN_CLANG_FORMAT_VERSION


def findClangFormat():
	"""Tries to find clang-format-XXX variants within the path"""
	cmd = "clang-format"
	allowed = [cmd] + [cmd + "-" + str(x) for x in range(nkt.MIN_CLANG_FORMAT_VERSION, nkt.MAX_CLANG_FORMAT_VERSION + 1)]
	for candidate in allowed:
		if isSupported(candidate):
			if nkt.isVerbose():
				version = str(subprocess.check_output([candidate, "--version"],
					universal_newlines=True, env=getEnv())).strip()
				print("clang-format: %s\n -> Version: %s" % (candidate, version))

			return candidate

	raise FileNotFoundError("clang-format binary not found. We searched for:\n " + "\n ".join(allowed))

def unsubscribedToFormat(filename, pattern = "no-networkit-format"):
	"""If pattern is present within the file, this file unsubscribed to auto formatting."""
	with open(filename, 'r') as file:
		return any( ((pattern in line) for line in file) )

def runClangFormat(inputFilename, outputFilename, clangFormatBinary = 'clang-format-8'):
	"""Execute clang-format onto inputFilename and stores the result in outputFilename"""
	with open(outputFilename, "w") as outfile:
		subprocess.call([clangFormatBinary, '-style=file', inputFilename], stdout=outfile)


nkt.setup()

# we need to set pwd in order for clang-format to find the .clang-format file!
os.chdir(nkt.getNetworKitRoot())

numberNonCompliant = 0
numberFileSkipped = 0

clangFormatCommand = findClangFormat()
with tempfile.TemporaryDirectory(dir=nkt.getNetworKitRoot()) as tempDir:
	files = nkt.getCXXFiles()
	for file in files:
		if unsubscribedToFormat(file):
			numberFileSkipped += 1
			continue

		tempFile = os.path.join(tempDir, 'cfOutput')
		runClangFormat(file, tempFile, clangFormatCommand)

		if not filecmp.cmp(file, tempFile, shallow=False):
			numberNonCompliant += 1
			nkt.reportChange(file + " is non-compliant")

			if nkt.doReportDiff():
				nkt.computeAndReportDiff(file, tempFile)

			if not nkt.isReadonly():
				os.replace(tempFile, file)

print("Scanned %d files (skipped %d files without subscription). Non-compliant files: %d." %
	  (len(files), numberFileSkipped, numberNonCompliant))

if numberNonCompliant > 0:
	nkt.failIfReadonly(__file__)
