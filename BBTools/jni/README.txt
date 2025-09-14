The accelerated versions of BBMap, Dedupe, BBMerge, and IceCreamFinder rely on the addition of a small amount C code that you may need to compile on your specific machine to take full advantage of its specific architecture.  However, compiled versions are already provided for Linux and OSX (libbbtoolsjni.so and libbbtoolsjni.dylib) with libbbtoolsjni.so compiled by Brian Bushnell on a Sandy Bridge machine and libbbtoolsjni.dylib compiled by Jie Wang on a Macbook; these seem to work fine on different architectures and operating systems.
To recompile, most C compilers should suffice. On Linux and OS X, we use gcc (gcc.gnu.org). The compiling process will create a library that Java can then load and use during execution. Simple makefiles for OSX and Linux have been provided. To compile the accelerated versions of BBTools on OS X or Linux, change your directory to "bbmap/jni" and refer to the following commands:

Linux:
make -f makefile.linux

OS X:
make -f makefile.osx

Windows:
If you are familiar with cmake and have it installed on your system, there is a CMakeLists.txt file that should get you most of the way there, but a full Windows build has not been tested at the moment.

After the "make" command, a "libbbtoolsjni.xx" file should appear in the "jni" directory. Once you have the libraries built, run BBMap, Dedupe, BBMerge, or IceCreamFinder with the addition of the "usejni" flag. If Java complains that it can't find the library you just compiled, the -Djava.library.path=<dir> flag in the bash scripts is what tells Java where to look for native library, and it should already be pointing to the "jni" directory. Java also looks for specific library suffixes on different operating systems, e.g. OSX: .dylib, Linux: .so, Windows: .dll.

Note:
Before building this, you must have a JDK installed (Java 8+) and the JAVA_HOME environment variable set (to something like /usr/bin/java/jdk-11.0.2/).
