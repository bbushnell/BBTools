package tax;

import java.io.File;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import parse.Parse;
import shared.Shared;
import shared.Tools;
import structures.IntList;

/**
 * Legacy one-dimensional GI-to-TaxID lookup table.
 * This predates {@link GiToTaxid}'s sparse two-dimensional arrays and keeps one
 * int slot for every GI from zero through the largest observed ID. Header helpers
 * can also pass explicit {@code ncbi|} or {@code tid|} TaxIDs through unchanged.
 * @author Brian Bushnell
 * @date Mar 10, 2015
 */
public class GiToTaxidInt{
	
	/** Constructs a legacy dense GI-to-TaxID utility instance. */
	public GiToTaxidInt(){}

	/** Builds and writes a serialized dense GI table.
	 * @param args Raw GI table path followed by the serialized output path */
	public static void main(String[] args){
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		ReadWrite.ZIPLEVEL=9;
		ReadWrite.PIGZ_BLOCKSIZE=256;
//		ReadWrite.PIGZ_ITERATIONS=30;
		
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			parse.Parser.parseZip(arg, a, b);
		}
//		if(args.length>2 && false){//Run a test
//			test(args);
//		}else 
		if(args.length>=2){//Write array
			initialize(args[0]);
			ReadWrite.write(array, args[1], true);
		}
	}
	
	/** Manual smoke test for an already-initialized GI table.
	 * @param args Optional taxonomy tree path in args[0], consulted only when
	 *        {@code args.length>1} */
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
	
	/** Parses a GI header and returns its stored dense-table value using {@code |}.
	 * @param s Header text, optionally beginning with {@code >}
	 * @return Stored TaxID or missing sentinel for an in-range GI, or -1 under
	 *         disabled assertions when the GI is absent or out of range */
	public static int parseGiToTaxid(String s){return parseGiToTaxid(s, '|');}

	/** Parses a GI header and returns its stored dense-table value.
	 * @param s Header text, optionally beginning with {@code >}
	 * @param delimiter Primary delimiter after {@code gi}
	 * @return Stored TaxID or missing sentinel for an in-range GI, or -1 under
	 *         disabled assertions when the GI is absent or out of range */
	public static int parseGiToTaxid(String s, char delimiter){
		int x=parseGiNumber(s, delimiter);
		assert(x>=0) : s;
		assert(array!=null) : "To use gi numbers, you must load a gi table.";
//		if(x>=array.length || array[x]<0){x=(int)(Math.random()*array.length);} //Test to make sure array is nonempty.
		if(x>=0 && x<array.length){return array[x];}
		assert(x<array.length) : "The GI number "+x+" is too big.\n"
				+ "Please update the gi table with the latest version from NCBI as per the instructions in gitable.sh.\n"
				+ "To ignore this problem, please run with the -da flag.\n";
		return -1;
	}
	

	/** Parses a byte-form GI header and returns its stored dense-table value using {@code |}.
	 * @param s Header bytes, optionally beginning with {@code >}
	 * @return Stored value for an in-range GI, or -1 for an absent GI header;
	 *         out-of-range GIs assert, or return -1 with assertions disabled */
	public static int parseGiToTaxid(byte[] s){return parseGiToTaxid(s, '|');}

	/** Parses a byte-form GI header and returns its stored dense-table value.
	 * @param s Header bytes, optionally beginning with {@code >}
	 * @param delimiter Primary delimiter after {@code gi}
	 * @return Stored value for an in-range GI, or -1 for an absent GI header;
	 *         out-of-range GIs assert, or return -1 with assertions disabled */
	public static int parseGiToTaxid(byte[] s, char delimiter){
		long x=parseGiNumber(s, delimiter);
		if(x>=0 && x<array.length){return array[(int)x];}
		if(x<0){return -1;}
		assert(false) : x;
		return -1;
	}
	
	/** Parses a GI number from string headers.
	 * @param s Header text, optionally beginning with {@code >}
	 * @param delimiter Primary delimiter after {@code gi}
	 * @return Parsed GI number, or -1 when no GI header is present */
	static int parseGiNumber(String s, char delimiter){
		if(s==null || s.length()<4){return -1;}
//		System.err.println("a");
		// Strip one leading '>' and keep parsing a GI number; do not call
		// getID(), which maps a parsed GI number through the table.
		if(s.charAt(0)=='>'){return parseGiNumber(s.substring(1), delimiter);}
//		System.err.println("b");
		if(!s.startsWith("gi")){return -1;}
//		System.err.println("c");
//		System.err.println("d");
		int initial=s.indexOf(delimiter);
//		System.err.println("e");
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
//			System.err.println("f");
//			System.err.println("g");
		}
//		System.err.println("h");
		if(!Tools.isDigit(s.charAt(initial+1))){return -1;}
//		System.err.println("i");
		
		int number=0;
		for(int i=initial+1; i<s.length(); i++){
			char c=s.charAt(i);
			if(c==delimiter){break;}
			assert(Tools.isDigit(c));
			number=(number*10)+(c-'0');
		}
//		System.err.println("j: "+number);
		return number;
	}
	
	/** Parses an explicit {@code ncbi} or {@code tid} TaxID header.
	 * @param s Header text, optionally beginning with {@code >}
	 * @param delimiter Primary delimiter after the prefix
	 * @return Parsed TaxID, or -1 for an unrecognized prefix, invalid first digit,
	 *         or value above {@link Integer#MAX_VALUE}. Later invalid digits also
	 *         return -1 when assertions are disabled.
	 * @throws AssertionError If assertions are enabled and a non-digit is
	 *         encountered after the first digit, before the field terminator */
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
	

	/** Returns the explicit TaxID or stored legacy GI value in a String header.
	 * @param s Header text, optionally beginning with {@code >}
	 * @return Explicit TaxID, stored TaxID or missing sentinel for an in-range GI,
	 *         or -1 when neither prefix is present */
	public static int getID(String s){return getID(s, '|');}

	/** Returns the explicit TaxID or stored legacy GI value in a String header.
	 * @param s Header text, optionally beginning with {@code >}
	 * @param delimiter Primary delimiter after the header prefix
	 * @return Explicit TaxID, stored TaxID or missing sentinel for an in-range GI,
	 *         or -1 when neither prefix is present */
	public static int getID(String s, char delimiter){
		int x=parseTaxidNumber(s, delimiter);
		if(x>=0){return x;}
		x=parseGiNumber(s, delimiter);
		//TODO: Confirm the intended range policy: this indexes the table directly,
		// unlike parseGiToTaxid(), which returns -1 for oversized GIs under -da.
		if(x>=0){return array[x];}
		return -1;
	}
	
	/** Parses a GI number from byte-form headers.
	 * @param s Header bytes, optionally beginning with {@code >}
	 * @param delimiter Primary delimiter after {@code gi}
	 * @return Parsed GI number, or -1 when no GI header is present */
	static int parseGiNumber(byte[] s, char delimiter){
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
		return (int)number;
	}
	
	/** Parses an explicit {@code ncbi} or {@code tid} TaxID from byte-form headers.
	 * @param s Header bytes, optionally beginning with {@code >}
	 * @param delimiter Primary delimiter after the prefix
	 * @return Parsed TaxID, or -1 for an unrecognized prefix, invalid first digit,
	 *         or value above {@link Integer#MAX_VALUE}. Later invalid digits also
	 *         return -1 when assertions are disabled.
	 * @throws AssertionError If assertions are enabled and a non-digit is
	 *         encountered after the first digit, before the field terminator */
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

	/** Returns the stored legacy GI value or explicit TaxID in a byte-form header.
	 * @param s Header bytes, optionally beginning with {@code >}
	 * @return Stored TaxID or missing sentinel for an in-range GI, explicit TaxID,
	 *         or -1 when neither prefix is present */
	public static int getID(byte[] s){return getID(s, '|');}

	/** Returns the stored legacy GI value or explicit TaxID in a byte-form header.
	 * @param s Header bytes, optionally beginning with {@code >}
	 * @param delimiter Primary delimiter after the header prefix
	 * @return Stored TaxID or missing sentinel for an in-range GI, explicit TaxID,
	 *         or -1 when neither prefix is present */
	public static int getID(byte[] s, char delimiter){
		int x=parseGiNumber(s, delimiter);
		//TODO: Confirm the intended range policy: this indexes the table directly,
		// unlike parseGiToTaxid(), which returns -1 for oversized GIs under -da.
		if(x>=0){return array[x];}
		return parseNcbiNumber(s, delimiter);
	}
	
	/**
	 * Looks up the stored value for one GI number by direct array index.
	 * @param gi GI number to look up
	 * @return Stored TaxID or missing sentinel from the dense table
	 */
	public static int getID(long gi){
		assert(gi>=0) : gi;
		assert(gi<Integer.MAX_VALUE) : gi+" > "+Integer.MAX_VALUE;
		assert(gi<array.length) : gi+", "+array.length;
		return array[(int)gi];
	}
	
	/** Loads a serialized one-dimensional array or builds one from GI input rows.
	 * Repeated calls with the same file name leave existing static state intact.
	 * @param fname Serialized .int1d table or raw GI table input */
	public static void initialize(String fname){
		assert(fname!=null);
		if(fileString==null || !fileString.equals(fname)){
			//TODO: Probable bug - this locks GiToTaxid.class while publishing
			// GiToTaxidInt's array, fileString and initialized state.
			synchronized(GiToTaxid.class){
				if(!initialized || fileString==null || !fileString.equals(fname)){
					final int[] table=(fname.contains(".int1d") ?
							ReadWrite.read(int[].class, fname, true) : makeArray(fname));
					array=table;
					fileString=fname;
				}
				initialized=true;
			}
		}
	}
	
	/** Returns true after the dense GI table has been loaded.
	 * @return True when a dense GI table is resident */
	public static boolean isInitialized(){return initialized;}
	
	/** Clears the resident array and loader state. */
	public static synchronized void unload(){
		array=null;
		fileString=null;
		initialized=false;
	}
	
	/** Builds a dense GI-indexed array from one or more raw mapping files.
	 * @param fnames Single path, comma list, or relative wildcard using '#'
	 * @return Dense GI-indexed TaxID array */
	private static int[] makeArray(String fnames){
		String[] split;
		if(new File(fnames).exists()){split=new String[] {fnames};}
		else if(fnames.indexOf(',')>=0){split=fnames.split(",");}
		else if(fnames.indexOf('#')>=0){
			assert(fnames.indexOf("/")<0) : "Note: Wildcard # only works for relative paths in present working directory.";
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
		
		IntList list=new IntList();
//		assert(max<Integer.MAX_VALUE) : "Overflow.";
//		int[] x=new int[(int)max+1];
//		Arrays.fill(x, -1);
		
		long total=0;
		for(String s : split){
			long count=addToList(s, list);
			total+=count;
		}
		return list.shrink().array;
	}
	
	/** Adds valid GI rows from one raw mapping file into a dense IntList.
	 * @param fname Input mapping file
	 * @param list Destination indexed by GI number
	 * @return Count of valid digit-ending rows seen before range filtering */
	private static long addToList(String fname, IntList list){
		boolean warned=false;
		//TODO: Probable bug - this ByteFile is not closed if line parsing or
		// IntList insertion fails before the loop reaches its explicit close.
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		long count=0, invalid=0;
		byte[] line=bf.nextLine();
		while(line!=null){
			if(line.length>0 && Tools.isDigit(line[line.length-1])){//Invalid lines will end with tab or na
				count++;
				int tab2=Tools.indexOfNth(line, '\t', 2);
				int tab3=Tools.indexOfNth(line, '\t', 1, tab2+1);
				assert(tab2>0 && (tab2<tab3) && tab3<line.length) : tab2+", "+tab3+", "+line.length;
				assert(tab2<line.length && line[tab2]=='\t') : tab2+", "+tab3+", '"+new String(line)+"'";
				assert(tab3<line.length && line[tab3]=='\t') : tab2+", "+tab3+", '"+new String(line)+"'";
//				assert(false) : tab2+", "+tab3+", '"+new String(line)+"'";
				int tid=Parse.parseInt(line, tab2+1, tab3);
				int gi=parseRawGi(line, tab3+1, line.length);
				if(gi<0){//A gi over 2.5b was observed May 3, 2021.
					invalid++;
				}else{
				assert(gi>=0) : "tid="+tid+", gi="+gi+", line=\n'"+new String(line)+"'";
				int old=list.get(gi);
				assert(old==0 || old==tid) : "Contradictory entries for gi "+gi+": "+old+" -> "+tid+"\n'"+new String(line)+"'\ntab2="+tab2+", tab3="+tab3;
				
				list.set(gi, tid);
				
				//assert(x[gi]==-1 || x[gi]==ncbi) : "Contradictory entries for gi "+gi+": "+x[gi]+" -> "+ncbi;
//				if(x[gi]!=-1 && x[gi]!=ncbi){
//					if(!warned){
//						System.err.println("***WARNING*** For file "+fname+":\n"+
//								("Contradictory entries for gi "+gi+": mapped to both taxID "+x[gi]+" and taxID "+ncbi)+
//								"\nThis may be an error from NCBI and you may wish to report it, but it is\n"
//								+ "being suppressed because NCBI data is known to contain multi-mapped gi numbers,\n"
//								+ "at least between nucleotide and protein, and gi numbers are deprecated anyway.");
//						warned=true;
//					}
//				}else{
//					x[gi]=ncbi;
//				}
				}
			}else{
				if(line.length==0){System.err.println(fname+", "+count);}//debug
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
	 * @return GI value, or -1 if the value exceeds the dense table limit */
	private static int parseRawGi(final byte[] line, final int start,
			final int stop){
		final int max=Shared.MAX_ARRAY_LEN-1;
		int number=0;
		int i=start;
		final int safeStop=i+Tools.min(stop-i, INT_SAFE_DIGITS);
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
	
//	private static int[] makeArrayOld(String fnames){
//		String[] split;
//		if(new File(fnames).exists()){split=new String[] {fnames};}
//		else{split=fnames.split(",");}
//		
//		long max=0;
//		for(String s : split){
//			max=Tools.max(max, findMaxID(s));
//		}
//		
//		assert(max<Integer.MAX_VALUE) : "Overflow.";
//		int[] x=new int[(int)max+1];
//		Arrays.fill(x, -1);
//		
//		long total=0;
//		for(String s : split){
//			long count=fillArray(s, x);
//			total+=count;
//		}
//		return x;
//	}
//	
//	private static long findMaxID(String fname){
//		ByteFile bf=ByteFile.makeByteFile(fname, true);
//		long count=0, max=0;
//		byte[] line=bf.nextLine();
//		while(line!=null){
//			count++;
//			int tab=Tools.indexOf(line, (byte)'\t');
//			long gi=Parse.parseLong(line, 0, tab);
//			max=Tools.max(max, gi);
//			line=bf.nextLine();
//		}
//		bf.close();
//		return max;
//	}
//	
//	private static long fillArray(String fname, int[] x){
//		boolean warned=false;
//		ByteFile bf=ByteFile.makeByteFile(fname, true);
//		long count=0;
//		byte[] line=bf.nextLine();
//		while(line!=null){
//			count++;
//			int tab=Tools.indexOf(line, (byte)'\t');
//			int gi=Parse.parseInt(line, 0, tab);
//			int ncbi=Parse.parseInt(line, tab+1, line.length);
//			//assert(x[gi]==-1 || x[gi]==ncbi) : "Contradictory entries for gi "+gi+": "+x[gi]+" -> "+ncbi;
//			if(x[gi]!=-1 && x[gi]!=ncbi){
//				if(!warned){
//					System.err.println("***WARNING*** For file "+fname+":\n"+
//							("Contradictory entries for gi "+gi+": mapped to both taxID "+x[gi]+" and taxID "+ncbi)+
//							"\nThis may be an error from NCBI and you may wish to report it, but it is\n"
//							+ "being suppressed because NCBI data is known to contain multi-mapped gi numbers,\n"
//							+ "at least between nucleotide and protein, and gi numbers are deprecated anyway.");
//					warned=true;
//				}
//			}else{
//				x[gi]=ncbi;
//			}
//			line=bf.nextLine();
//		}
//		if(verbose){System.err.println("Count: "+count);}
//		bf.close();
//		return count;
//	}
	
	/** Dense resident GI-indexed TaxID table. */
	private static int[] array;
	/** Path used to populate the resident table. */
	private static String fileString;
	
	/** Enables verbose load counters. */
	public static boolean verbose=false;
	/** True after a successful table load and before unload. */
	private static boolean initialized=false;
	/** Longest decimal prefix that fits in an int without overflow checks. */
	private static final int INT_SAFE_DIGITS=9;
}
