# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.Staurn.Debug:
/Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/Debug/Staurn:
	/bin/rm -f /Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/Debug/Staurn


PostBuild.Staurn.Release:
/Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/Release/Staurn:
	/bin/rm -f /Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/Release/Staurn


PostBuild.Staurn.MinSizeRel:
/Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/MinSizeRel/Staurn:
	/bin/rm -f /Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/MinSizeRel/Staurn


PostBuild.Staurn.RelWithDebInfo:
/Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/RelWithDebInfo/Staurn:
	/bin/rm -f /Users/Charles/Library/Mobile\ Documents/com~apple~CloudDocs/Project/Graphics_from_AtoZ/RelWithDebInfo/Staurn




# For each target create a dummy ruleso the target does not have to exist
