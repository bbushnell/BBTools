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
 * Maps GenInfo identifiers to NCBI taxonomy IDs using efficient 2D array storage.
 * Parses gi numbers from FASTA headers and converts them to taxonomy IDs for taxonomic
 * classification. Uses bit-shifting to divide gi space into upper/lower components
 * for memory-efficient sparse storage.
 *
 * @author Brian Bushnell
 * @date Mar 10, 2015
 */
public class GiToTaxid {
	
	public static void main(String[] args){
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		ReadWrite.ZIPLEVEL=9;
		ReadWrite.PIGZ_BLOCKSIZE=256;
//		ReadWrite.PIGZ_ITERATIONS=30;
		
		String names=args[0];
		String out=args.length>1 ? args[1] : null;
		for(int i=1; i<args.length; i++) {
			if(args[i].contains(".accession2taxid")) {
				names=names+","+args[i];
//			}else if(args[i].contains("assembly_summary_")) {
//				names=names+","+args[i];
			}else if(args[i].startsWith("gitable.")) {
				out=args[i];
			}
		}
		
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
			initialize(names);
			ReadWrite.write(array, out, true);
		}
	}
	
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
	
	public static int parseGiToTaxid(String s){return parseGiToTaxid(s, '|');}
	public static int parseGiToTaxid(String s, char delimiter){
		long x=parseGiNumber(s, delimiter);
		assert(x>=0) : x+", "+s;
		return getID(x);
	}
	

	public static int parseGiToTaxid(byte[] s){return parseGiToTaxid(s, '|');}
	public static int parseGiToTaxid(byte[] s, char delimiter){
		long x=parseGiNumber(s, delimiter);
		return x<0 ? -1 : getID(x);
	}
	
	static long parseGiNumber(String s, char delimiter){
		if(s==null || s.length()<4){return -1;}
		if(s.charAt(0)=='>'){return getID(s.substring(1), delimiter);}
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
		for(int i=initial+1; i<s.length(); i++){
			char c=s.charAt(i);
			if(c==delimiter || c==' '){break;}
			assert(Tools.isDigit(c)) : c+"\n"+s;
			number=(number*10)+(c-'0');
		}
		return number;
	}
	

	public static int getID(String s){return getID(s, '|');}
	public static int getID(String s, char delimiter){
		long x=parseTaxidNumber(s, delimiter);
		if(x>=0){return (int)x;}
		x=parseGiNumber(s, delimiter);
		return x<0 ? -1 : getID(x);
	}
	
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
		for(int i=initial+1; i<s.length; i++){
			byte c=s[i];
			if(c==delimiter){break;}
			assert(Tools.isDigit(c));
			number=(number*10)+(c-'0');
		}
		return number;
	}

	public static int getID(byte[] s){return getID(s, '|');}
	public static int getID(byte[] s, char delimiter){
		long x=parseGiNumber(s, delimiter);
		if(x>=0){return getID(x, true);}
		return parseNcbiNumber(s, delimiter);
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
	 * Looks up taxonomy ID for a gi number using 2D array storage.
	 * Uses bit-shifting to split gi into upper array index and lower position.
	 * Upper 34 bits select array slice, lower 30 bits select position within slice.
	 * @param gi GenInfo identifier to look up
	 * @param assertInRange Whether to assert on out-of-range gi numbers
	 * @return Taxonomy ID (0 if not present, -1 if invalid, -2 if out of range)
	 */
	public static int getID(long gi, boolean assertInRange){
		assert(initialized) : "To use gi numbers, you must load a gi table.";
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
	
	public static void initialize(String fnames){
		assert(fnames!=null);
		if(fileString==null || !fileString.equals(fnames)){
			synchronized(GiToTaxid.class){
				if(!initialized || fileString==null || !fileString.equals(fnames)){
					fileString=fnames;
					if(fnames.contains(".int2d")){
						array=ReadWrite.read(int[][].class, fnames, true);
						maxGiLoaded=-1;
						if(array!=null && array.length>0){
							int upper=array.length-1;
							int[] section=array[upper];
							int lower=section.length-1;
							maxGiLoaded=(((long)upper)<<SHIFT)|lower;
						}
					}else if(fnames.contains(".int1d")){
						throw new RuntimeException("Old gi table format filename "+fnames+".\n"
								+ "Current files should end in .int2d.");
						
					}else{
						array=makeArray(fnames);
					}
				}
				initialized=true;
			}
		}
	}
	
	public static boolean isInitialized(){return initialized;}
	
	public static synchronized void unload(){
		maxGiLoaded=-1;
		array=null;
		fileString=null;
		initialized=false;
	}
	
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
			table[i]=lists[i].array;
		}
		return table;
	}
	
	private static long addToList(String fname, IntList[] lists){
		boolean warned=false;
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		long count=0, invalid=0;
		byte[] line=bf.nextLine();
		while(line!=null){
			if(line.length>0 && Tools.isDigit(line[line.length-1])){//Invalid lines will end with tab or na
				count++;
				int tab2=-1;
				try{
					tab2=Tools.indexOfNth(line, '\t', 2);
				}catch(Throwable e){
					System.err.println(fname+", '"+new String(line)+"'");
					throw new RuntimeException(e);
				}
				int tab3=Tools.indexOfNth(line, '\t', 1, tab2+1);
				assert(tab2>0 && (tab2<tab3) && tab3<line.length) : tab2+", "+tab3+", "+line.length;
				assert(tab2<line.length && line[tab2]=='\t') : tab2+", "+tab3+", '"+new String(line)+"'";
				assert(tab3<line.length && line[tab3]=='\t') : tab2+", "+tab3+", '"+new String(line)+"'";
				//assert(false) : tab2+", "+tab3+", '"+new String(line)+"'";
				int tid=Parse.parseInt(line, tab2+1, tab3);
				int gi=Parse.parseInt(line, tab3+1, line.length);
				if(gi<0){
					invalid++;
				}else{
					assert(gi>=0) : "tid="+tid+", gi="+gi+", line=\n'"+new String(line)+"'";
					int old=setID(gi, tid, lists);
//					assert(old<1 || old==tid) : "Contradictory entries for gi "+gi+": "+old+" -> "+tid+"\n'"+new String(line)+"'\ntab2="+tab2+", tab3="+tab3;
					if(old>0 && old!=tid) {
						if(!warned) {
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
	
	private static int getID(long gi, IntList[] lists){
		assert(gi>=0) : "gi number "+gi+" is invalid.";
		final long upper=gi>>>SHIFT;
		final int lower=(int)(gi&LOWERMASK);
		assert(upper<Shared.MAX_ARRAY_LEN) : gi+", "+upper;
		IntList list=lists[(int)upper];
		return lower<0 ? -1 : lower>=list.size ? -2 : list.get(lower);
	}
	
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
	
	private static long contradictions=0;
	private static long maxGiLoaded=-1;
	private static int[][] array;
	private static final int SHIFT=30;
	private static final long UPPERMASK=(-1L)<<SHIFT;
	private static final long LOWERMASK=~UPPERMASK;
	
	private static String fileString;
	
	public static boolean verbose=false;
	private static boolean initialized=false;
}
