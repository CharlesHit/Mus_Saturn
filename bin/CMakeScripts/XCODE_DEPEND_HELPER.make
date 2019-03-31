# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.SATURN.Debug:
/Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/Debug/SATURN:
	/bin/rm -f /Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/Debug/SATURN


PostBuild.SATURN.Release:
/Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/Release/SATURN:
	/bin/rm -f /Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/Release/SATURN


PostBuild.SATURN.MinSizeRel:
/Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/MinSizeRel/SATURN:
	/bin/rm -f /Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/MinSizeRel/SATURN


PostBuild.SATURN.RelWithDebInfo:
/Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/RelWithDebInfo/SATURN:
	/bin/rm -f /Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/RelWithDebInfo/SATURN




# For each target create a dummy ruleso the target does not have to exist
