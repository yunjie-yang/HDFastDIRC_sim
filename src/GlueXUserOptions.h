//
//
// This is borrowed directly from HDGeant4/src/GlueXUserOptions.* minus the Geant4 part
//
//

#ifndef _GLUEXUSEROPTIONS_
#define _GLUEXUSEROPTIONS_

#include <stdlib.h>
#include <string>
#include <list>
#include <map>

class GlueXUserOptions
{
 // General class for reading options for controlling the HDGeant4 simulation
 // from input configuration files, including a geant3-style control.in file.
 // The class should be extended in the future to accept input from other
 // sources, eg. a database, by adding new ReadXXX methods.

 public:
   GlueXUserOptions();
   //GlueXUserOptions(const GlueXUserOptions &src);
   //GlueXUserOptions &operator=(const GlueXUserOptions &src);
   ~GlueXUserOptions();

   //static GlueXUserOptions *GetInstance();

   int ReadControl_in(const char *ctrlin="control.in");

   int Find(const char *name, std::map<int, std::string> &value) const;
   int Find(const char *name, std::map<int, double> &value) const;
   int Find(const char *name, std::map<int, int> &value) const;

 private:
   std::string askey(const std::string name) const;

   static std::list<GlueXUserOptions*> fInstance;
   std::map<std::string, std::string> fKeyValue;
};

#endif
