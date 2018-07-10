import numpy as np
import matplotlib.pyplot as plt

class Visualiser:
    def __init__(self):
        pass

    def hist(self, filaments, att, title="", labels="Data", saveas="_hist.png"):
        """
            filaments: Instance of class Filaments or list of it
            att: attribute to plot
        """
        plt.figure()
        if isinstance(filaments, list):
            if not isinstance(labels, list):
                print "Labels should be a list too!"
                print "This might give strange results."
                # TODO do a fix

            for l, f in zip(labels, filaments):
                # Attribute check
                if not (att in f.atts()):
                    print f, "does not contain the requested attribute."
                    print "Available attributes:"
                    print f.atts()
                    print "Skipping this filament for plotting."

                plt.hist(np.log10(f.data()[att]), bins=50, label=l, normed=True)

        else:
            if not (att in filaments.atts()):
                print "Filaments object does not contain the requested attribute."
                print "Available attributes:"
                print f.atts()
                plt.close()
                return

            plt.hist(np.log10(filaments.data()[att]), bins=50, label=labels, normed=True)

        plt.legend(loc="best")
        plt.minorticks_on()
        plt.title(title)
        plt.xlabel("log10" + att)
        plt.tight_layout()
        plt.savefig(att + saveas)
        plt.close()
