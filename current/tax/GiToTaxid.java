package tax;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import shared.Shared;
import shared.Tools;
import structures.IntList;

/**
 * Maps GenInfo identifiers to NCBI taxonomy IDs through either a sparse
 * two-dimensional array or a disk-backed table. FASTA helper methods can also
 * pass explicit {@code ncbi|} and {@code tid|} TaxIDs through unchanged.
 *
 * @author Brian Bushnell
 * @date Mar 10, 2015
 */
public class GiToTaxid{
	
	/** Constructs a GI-to-TaxID utility instance. */
	public GiToTaxid(){}

	/**
	 * With at least two arguments, builds and writes a serialized in-memory GI
	 * table from {@code args[0]}, preserving comma-list input syntax and later
	 * {@code .accession2taxid} paths. Recognized compression flags are applied,
	 * recognized Java wrapper flags are ignored, and one non-input positional
	 * token supplies the output path.
	 * @param args GI/accession source paths, compression flags, and output path
	 */
	public static void main(String[] args){
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		ReadWrite.ZIPLEVEL=9;
		ReadWrite.PIGZ_BLOCKSIZE=256;
//		ReadWrite.PIGZ_ITERATIONS=30;
		
		String names=args[0];
		String out=null;
		for(int i=1; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=", 2);
			final String a=split[0].toLowerCase();
			final String b=split.length>1 ? split[1] : null;
			if(Parser.parseZip(arg, a, b) || isRecognizedWrapperFlag(a)){
				//Recognized CLI flag.
			}else if(arg.contains(".accession2taxid")){
				names=names+","+arg;
			}else if(arg.startsWith("-")){
				throw new IllegalArgumentException("Unknown option: "+arg);
			}else if(out==null){
				out=arg;
			}else{
				throw new IllegalArgumentException("Multiple output paths: "+
						out+" and "+arg);
			}
		}
		if(args.length>=2){//Write array
			if(out==null){
				throw new IllegalArgumentException("Missing output path.");
			}
			checkInputOutputAliases(names, out);
			initialize(names);
			ReadWrite.write(array, out, true);
		}
	}

	/**
	 * Reports whether an argument was consumed by the shell wrapper before Java
	 * launched when the program is run normally.
	 * @param a Lowercase argument name
	 * @return True for recognized Java-wrapper-only flags
	 */
	private static boolean isRecognizedWrapperFlag(String a){
		return a.equals("-ea") || a.equals("ea") || a.equals("-da") ||
				a.equals("da") || a.equals("-eoom") || a.equals("eoom") ||
				a.equals("-exitonoutofmemoryerror") ||
				a.equals("exitonoutofmemoryerror") || a.startsWith("-xmx") ||
				a.startsWith("-xms");
	}

	/**
	 * Rejects a table output that would overwrite one of the raw inputs.
	 * Validation follows {@link #makeArray(String)} path expansion so literal
	 * comma filenames and local {@code #} wildcards are compared as the loader
	 * will read them.
	 * @param names GI/accession input expression
	 * @param out Serialized output path
	 */
	private static void checkInputOutputAliases(String names, String out){
		final File output=new File(out);
		final ArrayList<File> inputs=inputFiles(names);
		for(File input : inputs){
			if(sameFile(input, output)){
				throw new IllegalArgumentException(
						"Input and output files alias: "+input+" and "+out);
			}
		}
	}

	/**
	 * Resolves one GI input expression exactly as {@link #makeArray(String)}
	 * chooses raw table paths.
	 * @param fnames Single path, comma-delimited paths, or a local {@code #}
	 * wildcard
	 * @return Ordered input files that the raw loader will open
	 */
	private static ArrayList<File> inputFiles(String fnames){
		final String[] split;
		if(new File(fnames).exists()){split=new String[] {fnames};}
		else if(fnames.indexOf(',')>=0){split=fnames.split(",");}
		else if(fnames.indexOf('#')>=0 && fnames.indexOf('/')<0){
			final File dir=new File(System.getProperty("user.dir"));
			final String prefix=fnames.substring(0, fnames.indexOf('#'));
			final String suffix=fnames.substring(fnames.indexOf('#')+1);
			final File[] array=dir.listFiles();
			final ArrayList<File> files=new ArrayList<File>();
			if(array!=null){
				for(File f : array){
					final String s=f.getName();
					if(s.startsWith(prefix) && s.endsWith(suffix)){
						files.add(new File(s));
					}
				}
			}
			return files;
		}else{split=new String[] {fnames};}
		final ArrayList<File> files=new ArrayList<File>(split.length);
		for(String s : split){files.add(new File(s));}
		return files;
	}

	/**
	 * Tests two paths for identity through canonical names and same-file
	 * metadata when both paths already exist.
	 * @param a First path
	 * @param b Second path
	 * @return True if the paths identify the same file
	 */
	private static boolean sameFile(File a, File b){
		final File ca=canonicalFile(a);
		final File cb=canonicalFile(b);
		if(ca.equals(cb)){return true;}
		if(ca.exists() && cb.exists()){
			try{
				return Files.isSameFile(ca.toPath(), cb.toPath());
			}catch(IOException e){
				throw new RuntimeException(
						"Failed to compare GI table path identity: "+a+
						" and "+b, e);
			}
		}
		return false;
	}

	/**
	 * Canonicalizes a path for alias checks.
	 * @param path Path to canonicalize
	 * @return Canonical file
	 */
	private static File canonicalFile(String path){
		return canonicalFile(new File(path));
	}

	/**
	 * Canonicalizes a file for alias checks.
	 * @param file File to canonicalize
	 * @return Canonical file
	 */
	private static File canonicalFile(File file){
		try{
			return file.getCanonicalFile();
		}catch(IOException e){
			throw new RuntimeException(
					"Failed to canonicalize GI table path: "+file, e);
		}
	}
	
	/**
	 * Manual smoke test for an already initialized GI table and optional tree.
	 * Only {@code args[0]} is consumed, and only when {@code args.length>1};
	 * this helper never loads GI table arguments itself.
	 * @param args Optional taxonomy tree path in position 0
	 */
	public static void test(String[] args){
		System.err.println(getID(1000));
		System.err.println(getID(10000));
		System.err.println(getID(10001));
		System.err.println(getID(10002));
		System.err.println(getID(10003));
		System.err.println(getID(10004));
		System.err.println(getID(10005));
		System.err.println(getID(100000));
		System.err.println(getID(1000000));
		System.err.println(getID(10000000));
		
		TaxTree tree=null;
		if(args.length>1){
			tree=TaxTree.loadTaxTree(args[0], System.err, true, true);
		}
		
		System.err.println("Strings:");
		int x;
		x=getID("gi|18104025|emb|AJ427095.1| Ceratitis capitata centromeric or pericentromeric satellite DNA, clone 44");
		System.err.println(x);
		if(tree!=null){
			System.err.println(tree.getNode(x));
			tree.incrementRaw(x, 30);
		}
		x=getID("gi|15982920|gb|AY057568.1| Arabidopsis thaliana AT5g43500/MWF20_22 mRNA, complete cds");
		System.err.println(x);
		if(tree!=null){
			System.err.println(tree.getNode(x));
			tree.incrementRaw(x, 40);
		}
		x=getID("gi|481043749|gb|KC494054.1| Plesiochorus cymbiformis isolate ST05-58 internal transcribed spacer 2, partial sequence");
		System.err.println(x);
		if(tree!=null){
			System.err.println(tree.getNode(x));
			tree.incrementRaw(x, 20);
		}
		
		if(tree!=null){
			tree.percolateUp();
			ArrayList<TaxNode> nodes=tree.gatherNodesAtLeastLimit(35);
			for(TaxNode n : nodes){
				System.err.println(n);
			}
		}
	}
	
	/**
	 * Parses a pipe-delimited GI header and resolves it in the active table.
	 * @param s Header or GI token
	 * @return TaxID, or a negative lookup sentinel
	 */
	public static int parseGiToTaxid(String s){return parseGiToTaxid(s, '|');}
	/**
	 * Parses a GI header and resolves it in the active table.
	 * @param s Header or GI token
	 * @param delimiter Primary delimiter between the {@code gi} prefix and number
	 * @return TaxID, or a negative lookup sentinel
	 */
	public static int parseGiToTaxid(String s, char delimiter){
		long x=parseGiNumber(s, delimiter);
		assert(x>=0) : x+", "+s;
		return getID(x);
	}
	
	/**
	 * Parses a pipe-delimited GI byte header and resolves it in the active table.
	 * @param s Header or GI token bytes
	 * @return TaxID, or a negative lookup sentinel
	 */
	public static int parseGiToTaxid(byte[] s){return parseGiToTaxid(s, '|');}
	/**
	 * Parses a GI byte header and resolves it in the active table.
	 * @param s Header or GI token bytes
	 * @param delimiter Primary delimiter between the {@code gi} prefix and number
	 * @return TaxID, or a negative lookup sentinel
	 */
	public static int parseGiToTaxid(byte[] s, char delimiter){
		long x=parseGiNumber(s, delimiter);
		return x<0 ? -1 : getID(x);
	}
	
	/**
	 * Extracts a GI number from a string token with fallback separators.
	 * @param s Header or GI token
	 * @param delimiter Primary separator to try before {@code ~} and {@code _}
	 * @return Parsed GI number, or -1 if the token is not a GI header
	 */
	static long parseGiNumber(String s, char delimiter){
		if(s==null || s.length()<4){return -1;}
		if(s.charAt(0)=='>'){return parseGiNumber(s.substring(1), delimiter);}
		if(!s.startsWith("gi")){return -1;}
		int initial=s.indexOf(delimiter);
		if(initial<0){
			if(delimiter!='~'){
				delimiter='~';
				initial=s.indexOf(delimiter);
			}
			if(initial<0){
				delimiter='_';
				initial=s.indexOf(delimiter);
			}
			if(initial<0){return -1;}
		}
		if(!Tools.isDigit(s.charAt(initial+1))){return -1;}
		
		long number=0;
		for(int i=initial+1; i<s.length(); i++){
			char c=s.charAt(i);
			if(c==delimiter){break;}
			assert(Tools.isDigit(c));
			number=(number*10)+(c-'0');
		}
		return number;
	}
	
	/**
	 * Extracts an NCBI taxonomy ID from an {@code ncbi} or {@code tid} string.
	 * @param s Header or taxID token
	 * @param delimiter Primary separator before falling back to {@code _}
	 * @return Parsed taxonomy ID, or -1 for an unrecognized prefix, invalid first
	 * digit, or value above {@link Integer#MAX_VALUE}. Later invalid digits also
	 * return -1 when assertions are disabled.
	 * @throws AssertionError If assertions are enabled and a non-digit is
	 * encountered after the first digit, before the field terminator
	 */
	public static int parseTaxidNumber(String s, char delimiter){
		if(s==null || s.length()<5){return -1;}
		if(s.charAt(0)=='>'){return parseTaxidNumber(s.substring(1), delimiter);}
		if(!s.startsWith("ncbi") && !s.startsWith("tid")){return -1;}
		int initial=s.indexOf(delimiter);
		if(initial<0){
			delimiter='_';
			initial=s.indexOf(delimiter);
			if(initial<0){return -1;}
		}
		if(!Tools.isDigit(s.charAt(initial+1))){return -1;}
		
		int number=0;
		int i=initial+1;
		final int safeStop=i+Tools.min(s.length()-i, INT_SAFE_DIGITS);
		for(; i<safeStop; i++){
			char c=s.charAt(i);
			if(c==delimiter || c==' '){return number;}
			assert(Tools.isDigit(c)) : c+"\n"+s;
			if(!Tools.isDigit(c)){return -1;}
			number=(number*10)+(c-'0');
		}
		return parseTaxidNumberTail(s, delimiter, i, number);
	}

	/** Continues explicit String TaxID parsing after the first safe int digits.
	 * @param s Header text already accepted by {@link #parseTaxidNumber(String, char)}
	 * @param delimiter Active delimiter, including fallback to {@code _}
	 * @param i First unparsed digit offset
	 * @param number Value accumulated before this tail
	 * @return Parsed TaxID, or -1 for malformed digits under disabled assertions
	 *         or values above {@link Integer#MAX_VALUE}
	 * @throws AssertionError If assertions are enabled and a non-digit is
	 *         encountered before the field terminator */
	private static int parseTaxidNumberTail(String s, char delimiter,
			int i, int number){
		long longNumber=number;
		for(; i<s.length(); i++){
			char c=s.charAt(i);
			if(c==delimiter || c==' '){break;}
			assert(Tools.isDigit(c)) : c+"\n"+s;
			if(!Tools.isDigit(c)){return -1;}
			longNumber=(longNumber*10)+(c-'0');
			if(longNumber>Integer.MAX_VALUE){return -1;}
		}
		return longNumber>Integer.MAX_VALUE ? -1 : (int)longNumber;
	}
	
	/**
	 * Resolves a pipe-delimited string header as an explicit taxID or GI number.
	 * @param s Header or identifier token
	 * @return TaxID, or a negative lookup sentinel
	 */
	public static int getID(String s){return getID(s, '|');}
	/**
	 * Resolves a string header as an explicit taxID or GI number.
	 * @param s Header or identifier token
	 * @param delimiter Primary separator
	 * @return TaxID, or a negative lookup sentinel
	 */
	public static int getID(String s, char delimiter){
		long x=parseTaxidNumber(s, delimiter);
		if(x>=0){return (int)x;}
		x=parseGiNumber(s, delimiter);
		return x<0 ? -1 : getID(x);
	}
	
	/**
	 * Extracts a GI number from byte header data.
	 * @param s Header or GI token bytes
	 * @param delimiter Primary separator before falling back to {@code _}
	 * @return Parsed GI number, or -1 if the token is not a GI header
	 */
	static long parseGiNumber(byte[] s, char delimiter){
		if(s==null || s.length<4){return -1;}
		if(!Tools.startsWith(s, "gi") && !Tools.startsWith(s, ">gi")){return -1;}
		int initial=Tools.indexOf(s, (byte)delimiter);
		if(initial<0){
			delimiter='_';
			initial=Tools.indexOf(s, (byte)delimiter);
			if(initial<0){return -1;}
		}
		if(!Tools.isDigit(s[initial+1])){return -1;}
		
		long number=0;
		for(int i=initial+1; i<s.length; i++){
			byte c=s[i];
			if(c==delimiter){break;}
			assert(Tools.isDigit(c));
			number=(number*10)+(c-'0');
		}
		return number;
	}
	
	/**
	 * Extracts an NCBI taxonomy ID from byte data with an {@code ncbi} or
	 * {@code tid} prefix.
	 * @param s Header or taxID token bytes
	 * @param delimiter Primary separator before falling back to {@code _}
	 * @return Parsed taxonomy ID, or -1 for an unrecognized prefix, invalid first
	 * digit, or value above {@link Integer#MAX_VALUE}. Later invalid digits also
	 * return -1 when assertions are disabled.
	 * @throws AssertionError If assertions are enabled and a non-digit is
	 * encountered after the first digit, before the field terminator
	 */
	static int parseNcbiNumber(byte[] s, char delimiter){
		if(s==null || s.length<3){return -1;}
		if(!Tools.startsWith(s, "ncbi") && !Tools.startsWith(s, ">ncbi") && !Tools.startsWith(s, "tid") && !Tools.startsWith(s, ">tid")){return -1;}
		int initial=Tools.indexOf(s, (byte)delimiter);
		if(initial<0){
			delimiter='_';
			initial=Tools.indexOf(s, (byte)delimiter);
			if(initial<0){return -1;}
		}
		if(!Tools.isDigit(s[initial+1])){return -1;}
		
		int number=0;
		int i=initial+1;
		final int safeStop=i+Tools.min(s.length-i, INT_SAFE_DIGITS);
		for(; i<safeStop; i++){
			byte c=s[i];
			if(c==delimiter){return number;}
			assert(Tools.isDigit(c));
			if(!Tools.isDigit(c)){return -1;}
			number=(number*10)+(c-'0');
		}
		long longNumber=number;
		for(; i<s.length; i++){
			byte c=s[i];
			if(c==delimiter){break;}
			assert(Tools.isDigit(c));
			if(!Tools.isDigit(c)){return -1;}
			longNumber=(longNumber*10)+(c-'0');
			if(longNumber>Integer.MAX_VALUE){return -1;}
		}
		return longNumber>Integer.MAX_VALUE ? -1 : (int)longNumber;
	}

	/**
	 * Resolves a pipe-delimited byte header as a GI number or explicit taxID.
	 * @param s Header or identifier token bytes
	 * @return TaxID, or a negative lookup sentinel
	 */
	public static int getID(byte[] s){return getID(s, '|');}
	/**
	 * Resolves a byte header as a GI number or explicit taxID.
	 * @param s Header or identifier token bytes
	 * @param delimiter Primary separator
	 * @return TaxID, or a negative lookup sentinel
	 */
	public static int getID(byte[] s, char delimiter){
		long x=parseGiNumber(s, delimiter);
		if(x>=0){return getID(x, true);}
		return parseNcbiNumber(s, delimiter);
	}
	
	/**
	 * Loads disk-backed GI table for low-memory server mode.
	 * Callers must quiesce readers before replacing the active GI backend.
	 * @param path Path to the binary GI table file
	 */
	public static void initializeDisk(String path){
		DiskGiTable table=null;
		Throwable failure=null;
		try{
			table=new DiskGiTable(path);
			synchronized(GiToTaxid.class){
				final RuntimeException closeFailure=closeDiskTable();
				array=null;
				maxGiLoaded=-1;
				fileString=path;
				diskTable=table;
				DISK_MODE=true;
				initialized=true;
				table=null;
				if(closeFailure!=null){
					throw oldTableCloseFailure("publishing disk GI table", closeFailure);
				}
			}
			System.err.println("Loaded disk GI table: "+path);
		}catch(IOException e){
			final RuntimeException wrapped=new RuntimeException(
					"Failed to load disk GI table: "+path, e);
			failure=wrapped;
			throw wrapped;
		}catch(RuntimeException e){
			failure=e;
			throw e;
		}catch(Error e){
			failure=e;
			throw e;
		}finally{
			final RuntimeException closeFailure=close(table);
			if(closeFailure!=null){
				if(failure==null){throw closeFailure;}
				failure.addSuppressed(closeFailure);
			}
		}
	}

	/**
	 * Looks up taxonomy ID for a gi number in the loaded mapping table.
	 * @param gi GenInfo identifier to look up
	 * @return Taxonomy ID; 0 if not present, -1 if invalid (negative), -2 if out of range
	 */
	public static int getID(long gi){
		return getID(gi, true);
	}

	/**
	 * Looks up taxonomy ID for a gi number using 2D array storage or disk table.
	 * Uses bit-shifting to split gi into upper array index and lower position.
	 * Upper 34 bits select array slice, lower 30 bits select position within slice.
	 * @param gi GenInfo identifier to look up
	 * @param assertInRange Whether to assert on out-of-range gi numbers
	 * @return Taxonomy ID (0 if not present, -1 if invalid, -2 if out of range)
	 */
	public static int getID(long gi, boolean assertInRange){
		assert(initialized) : "To use gi numbers, you must load a gi table.";
		if(DISK_MODE){return diskTable.get(gi);}
		if(gi<0 || gi>maxGiLoaded){
			assert(!assertInRange) : gi<0 ? "gi number "+gi+" is invalid." :
				"The gi number "+gi+" is too big: Max loaded gi number is "+maxGiLoaded+".\n"
				+ "Please update the gi table with the latest version from NCBI"
				+ " as per the instructions in gitable.sh.\n"
				+ "To ignore this problem, please run with the -da flag.\n";
			return gi<0 ? -1 : -2;
		}
		final long upper=gi>>>SHIFT;
		final int lower=(int)(gi&LOWERMASK);
		assert(upper<Shared.MAX_ARRAY_LEN && upper<array.length) : gi+", "+upper+", "+array.length;
		final int[] slice=array[(int)upper];
		return slice==null || slice.length<=lower ? 0 : slice[lower];
	}
	
	/**
	 * Loads an in-memory GI table and makes it the active GI backend.
	 * Callers must quiesce readers before replacing the active GI backend.
	 * @param fnames Comma-separated or wildcard GI table input paths
	 */
	public static void initialize(String fnames){
		assert(fnames!=null);
		if(DISK_MODE || !initialized || fileString==null || !fileString.equals(fnames)){
			synchronized(GiToTaxid.class){
				if(DISK_MODE || !initialized || fileString==null || !fileString.equals(fnames)){
					int[][] table=null;
					long maxGi=-1;
					if(fnames.contains(".int2d")){
						table=ReadWrite.read(int[][].class, fnames, true);
						if(table!=null && table.length>0){
							int upper=table.length-1;
							int[] section=table[upper];
							int lower=section.length-1;
							maxGi=(((long)upper)<<SHIFT)|lower;
						}
					}else if(fnames.contains(".int1d")){
						throw new RuntimeException("Old gi table format filename "+fnames+".\n"
								+ "Current files should end in .int2d.");
						
					}else{
						final long oldMaxGi=maxGiLoaded;
						try{
							maxGiLoaded=-1;
							table=makeArray(fnames);
							maxGi=maxGiLoaded;
						}finally{
							maxGiLoaded=oldMaxGi;
						}
					}
					final RuntimeException closeFailure=closeDiskTable();
					array=table;
					maxGiLoaded=maxGi;
					fileString=fnames;
					DISK_MODE=false;
					initialized=true;
					if(closeFailure!=null){
						throw oldTableCloseFailure("publishing array GI table", closeFailure);
					}
				}
				initialized=true;
			}
		}
	}
	
	/**
	 * Reports whether an in-memory or disk GI backend is active.
	 * @return True after successful initialization and before unload
	 */
	public static boolean isInitialized(){return initialized;}
	
	/**
	 * Clears the active GI backend and closes a disk table if one is loaded.
	 * Callers must quiesce readers before unloading the active GI backend.
	 */
	public static synchronized void unload(){
		final RuntimeException closeFailure=closeDiskTable();
		maxGiLoaded=-1;
		array=null;
		fileString=null;
		DISK_MODE=false;
		initialized=false;
		if(closeFailure!=null){throw closeFailure;}
	}
	
	/**
	 * Builds an in-memory sparse GI table from one or more source files.
	 * @param fnames Single path, comma-delimited paths, or a local {@code #} wildcard
	 * @return Two-dimensional table indexed by upper/lower GI bits
	 */
	private static int[][] makeArray(String fnames){
		String[] split;
		if(new File(fnames).exists()){split=new String[] {fnames};}
		else if(fnames.indexOf(',')>=0){split=fnames.split(",");}
		else if(fnames.indexOf('#')>=0){
			assert(fnames.indexOf("/")<0) : "Note: Wildcard # only works for "
					+ "relative paths in present working directory.";
			File dir=new File(System.getProperty("user.dir"));
			String prefix=fnames.substring(0, fnames.indexOf('#'));
			String suffix=fnames.substring(fnames.indexOf('#')+1);
			
			File[] array=dir.listFiles();
			StringBuilder sb=new StringBuilder();
			String comma="";
			for(File f : array){
				String s=f.getName();
				if(s.startsWith(prefix) && s.endsWith(suffix)){
					sb.append(comma);
					sb.append(s);
					comma=",";
				}
			}
			split=sb.toString().split(",");
		}else{
			throw new RuntimeException("Invalid file: "+fnames);
		}
		
		int numLists=32;
		IntList[] lists=new IntList[numLists];
		
		long total=0;
		for(String s : split){
			long count=addToList(s, lists);
			total+=count;
		}
		System.err.println("total:          "+total);
		System.err.println("contradictions: "+contradictions);
		for(int i=0; i<lists.length; i++){
			if(lists[i]!=null && lists[i].size>0){
				lists[i].shrink();
				numLists=i+1;
			}
		}
		int[][] table=new int[numLists][];
		for(int i=0; i<numLists; i++){
			table[i]=(lists[i]==null ? new int[0] : lists[i].array);
		}
		return table;
	}
	
	/**
	 * Parses one GI source file into expandable upper-bit slices.
	 * @param fname Source GI mapping file
	 * @param lists Sparse upper-bit slices to fill
	 * @return Count of nonempty rows ending in a digit, before invalid-GI filtering
	 */
	private static long addToList(String fname, IntList[] lists){
		boolean warned=false;
		//TODO: Probable bug - this ignores ByteFile.close() status, and if
		// parsing or table insertion throws, the file never reaches close.
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		long count=0, invalid=0;
		byte[] line=bf.nextLine();
		while(line!=null){
			if(line.length>0 && Tools.isDigit(line[line.length-1])){//Invalid lines will end with tab or na
				count++;
				int tab2=-1, tab3=-1;
				try{
					tab2=Tools.indexOfNth(line, '\t', 2);
					tab3=Tools.indexOfNth(line, '\t', 1, tab2+1);
					if(tab3>0){
						assert(tab2>0 && (tab2<tab3) && tab3<line.length) : tab2+", "+tab3+", "+line.length;
						assert(tab2<line.length && line[tab2]=='\t') : tab2+", "+tab3+", '"+new String(line)+"'";
						assert(tab3<line.length && line[tab3]=='\t') : tab2+", "+tab3+", '"+new String(line)+"'";
					}
				}catch(Throwable e){
					System.err.println(fname+", '"+new String(line)+"'");
					throw new RuntimeException(e);
				}
				//assert(false) : tab2+", "+tab3+", '"+new String(line)+"'";
				int tid=Parse.parseInt(line, tab2+1, tab3>0 ? tab3 : line.length);
				long gi=(tab3>0 ? parseRawGi(line, tab3+1, line.length, lists) : -1);
				if(gi<0){
					invalid++;
				}else{
					assert(gi>=0) : "tid="+tid+", gi="+gi+", line=\n'"+new String(line)+"'";
					int old=setID(gi, tid, lists);
//					assert(old<1 || old==tid) : "Contradictory entries for gi "+gi+": "+old+" -> "+tid+"\n'"+new String(line)+"'\ntab2="+tab2+", tab3="+tab3;
					if(old>0 && old!=tid){
						if(!warned){
							System.err.println("Contradictory entries for gi "+gi+": "+old+" -> "+tid+"\n'"+new String(line)+"'\ntab2="+tab2+", tab3="+tab3);
						}
						warned=true;
						contradictions++;
					}
				}
			}else{
				//if(line.length==0){System.err.println(fname+", "+count);}//debug
				invalid++;
			}
			line=bf.nextLine();
		}
		if(verbose){System.err.println("Count: "+count+"; \tInvalid: "+invalid);}
		bf.close();
		return count;
	}

	/** Parses one raw-table GI field without allowing long overflow to wrap.
	 * @param line Complete raw input row
	 * @param start First GI digit offset
	 * @param stop Offset after the final GI digit
	 * @param lists Sparse upper-bit slices defining the supported raw GI range
	 * @return GI value, or -1 if the value exceeds the sparse table capacity */
	private static long parseRawGi(final byte[] line, final int start,
			final int stop, final IntList[] lists){
		final long max=(((long)lists.length)<<SHIFT)-1;
		long number=0;
		int i=start;
		final int safeStop=i+Tools.min(stop-i, LONG_SAFE_DIGITS);
		for(; i<safeStop; i++){
			final byte c=line[i];
			if(!Tools.isDigit(c)){return -1;}
			number=(number*10)+(c-'0');
		}
		for(; i<stop; i++){
			final byte c=line[i];
			if(!Tools.isDigit(c)){return -1;}
			final int digit=c-'0';
			if(number>(max-digit)/10){return -1;}
			number=(number*10)+digit;
		}
		return number>max ? -1 : number;
	}
	
	/**
	 * Looks up a GI while a new in-memory table is still under construction.
	 * @param gi GenInfo identifier
	 * @param lists Sparse upper-bit slices
	 * @return TaxID, -1 for invalid lower bits, or -2 past the loaded slice end
	 */
	private static int getID(long gi, IntList[] lists){
		assert(gi>=0) : "gi number "+gi+" is invalid.";
		final long upper=gi>>>SHIFT;
		final int lower=(int)(gi&LOWERMASK);
		assert(upper<Shared.MAX_ARRAY_LEN) : gi+", "+upper;
		IntList list=lists[(int)upper];
		return lower<0 ? -1 : lower>=list.size ? -2 : list.get(lower);
	}
	
	/**
	 * Stores one GI to taxID mapping in the expandable build table.
	 * @param gi GenInfo identifier
	 * @param tid NCBI taxonomy ID
	 * @param lists Sparse upper-bit slices
	 * @return Previous value at the selected GI, or an absent/out-of-range sentinel
	 */
	private static int setID(long gi, int tid, IntList[] lists){
		assert(gi>=0) : "gi number "+gi+" is invalid.";
		final long upper=gi>>>SHIFT;
		final int lower=(int)(gi&LOWERMASK);
		assert(upper<Shared.MAX_ARRAY_LEN) : gi+", "+upper;
		IntList list=lists[(int)upper];
		if(list==null){list=lists[(int)upper]=new IntList();}
		int old=lower<0 ? -1 : lower>=list.size ? -2 : list.get(lower);
		list.set(lower, tid);
		maxGiLoaded=Tools.max(gi, maxGiLoaded);
		return old;
	}

	/**
	 * Clears and closes the active disk table, returning cleanup failure to report
	 * after replacement state has been published.
	 * Readers must already be quiescent so no lookup can touch the old table.
	 * @return Close failure, or null after a clean close or no-op
	 */
	private static RuntimeException closeDiskTable(){
		final DiskGiTable table=diskTable;
		diskTable=null;
		return close(table);
	}

	/**
	 * Closes one disk table without mutating publication state.
	 * @param table Table to close; null is accepted
	 * @return Close failure, or null after a clean close or no-op
	 */
	private static RuntimeException close(final DiskGiTable table){
		if(table==null){return null;}
		try{
			table.close();
			return null;
		}catch(IOException e){
			return new RuntimeException("Failed to close disk GI table.", e);
		}
	}

	/**
	 * Wraps an old-table close failure with the replacement state already visible.
	 * @param action State transition already published
	 * @param closeFailure Failure returned by the old table close attempt
	 * @return Exception reporting both the close failure and current state
	 */
	private static RuntimeException oldTableCloseFailure(final String action,
			final RuntimeException closeFailure){
		return new RuntimeException("Failed to close previous disk GI table after "+
				action+"; initialized="+initialized+", DISK_MODE="+DISK_MODE+
				", fileString="+fileString+".", closeFailure);
	}
	
	/** True when GI lookups route to {@link #diskTable} instead of {@link #array}. */
	public static boolean DISK_MODE=false;
	/** Active disk-backed GI table, or null in in-memory mode. */
	private static DiskGiTable diskTable=null;
	/** Count of conflicting GI rows observed while building an in-memory table. */
	private static long contradictions=0;
	/** Largest GI present in the active or currently building in-memory table. */
	static long maxGiLoaded=-1;
	/** Sparse in-memory GI-to-taxID table indexed by upper/lower GI bits. */
	static int[][] array;
	/** Number of low GI bits used as an in-memory slice offset. */
	private static final int SHIFT=30;
	/** Mask selecting the upper slice-index bits of a GI. */
	private static final long UPPERMASK=(-1L)<<SHIFT;
	/** Mask selecting the lower offset bits of a GI. */
	private static final long LOWERMASK=~UPPERMASK;
	
	/** Path or path expression for the currently initialized GI backend. */
	private static String fileString;
	
	/** Enables verbose GI table-loading statistics. */
	public static boolean verbose=false;
	/** True after a GI backend has been loaded and published. */
	private static boolean initialized=false;
	/** Longest decimal prefix that fits in an int without overflow checks. */
	private static final int INT_SAFE_DIGITS=9;
	/** Longest decimal prefix that fits in a long without overflow checks. */
	private static final int LONG_SAFE_DIGITS=18;
}
