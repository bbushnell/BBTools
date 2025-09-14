package dna;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

import fileIO.ByteFile1;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Uses a ByteFile instead of TextFile for better speed and lower memory use.
 * @author Brian Bushnell
 * @date Jul 30, 2014
 *
 */
public class FastaToChromArrays2 {
	
//	Example:
//	jgi.FastaToChromArrays ecoli_K12.fa 1 writeinthread=false genscaffoldinfo=true retain waitforwriting=false
//	gzip=true chromc=false maxlen=536670912 writechroms=true minscaf=1 midpad=300 startpad=8000 stoppad=8000 nodisk=false
	
	public static void main(String[] args){
		main2(args);
	}
	
	public static ArrayList<ChromosomeArray> main2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		boolean oldWIT=WRITE_IN_THREAD;
		WRITE_IN_THREAD=true;
		
		String name=null;
		int genome=-1;
		int chroms=-1;
		String infile=null;
		boolean writeinfo=false;
		boolean genScaffoldInfo=true;
		boolean writeChroms=true;
		boolean scafprefixes=Data.scaffoldPrefixes;
		
		for(int i=0; i<args.length; i++){

			if(true){
				final String arg=args[i];
				final String[] split=arg.split("=");
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				
				if(a.equals("null")){
					//do nothing
				}else if(a.equals("path") || a.equals("root") || a.equals("tempdir")){
					Data.setPath(b);
				}else if(a.equals("name") || a.equals("organism")){
					name=b;
				}else if(a.equals("in") || a.equals("input") || a.equals("ref") || a.equals("fasta")){
					if(split.length<1 || "null".equalsIgnoreCase(b)){b=null;}
					infile=b;
				}else if(a.equals("build") || a.equals("genome")){
					genome=Integer.parseInt(b);
				}else if(a.equals("chroms")){
					chroms=Integer.parseInt(b);
				}else if(a.equals("writeinthread")){
					WRITE_IN_THREAD=Parse.parseBoolean(b);
				}else if(a.equals("nodisk")){
					NODISK=Parse.parseBoolean(b);
				}else if(a.equals("writeinfo")){
					writeinfo=Parse.parseBoolean(b);
				}else if(a.equals("padstart") || a.startsWith("startpad") || a.equals("padfront") || a.startsWith("frontpad")){
					START_PADDING=Integer.parseInt(b);
				}else if(a.equals("padstop") || a.startsWith("stoppad") || a.equals("padend") || a.startsWith("endpad")){
					END_PADDING=Integer.parseInt(b);
				}else if(a.equals("pad") || a.equals("padding")){
					START_PADDING=END_PADDING=Integer.parseInt(b);
				}else if(a.equals("midpad") || a.startsWith("padmid")){
					MID_PADDING=Integer.parseInt(b);
				}else if(a.startsWith("minscaf") || a.startsWith("mincontig")){
					MIN_SCAFFOLD=Integer.parseInt(b);
				}else if(a.equals("genscaffoldinfo")){
					genScaffoldInfo=Parse.parseBoolean(b);
					System.err.println("Set genScaffoldInfo="+genScaffoldInfo);
				}else if(a.equals("append") || a.equals("app")){
					append=Parse.parseBoolean(b);
				}else if(a.equals("overwrite") || a.equals("ow")){
					overwrite=Parse.parseBoolean(b);
				}else if(a.equals("mergescaffolds") || a.equals("mergecontigs") || (a.equals("merge"))){
					MERGE_SCAFFOLDS=Parse.parseBoolean(b);
					System.err.println("Set MERGE_SCAFFOLDS="+MERGE_SCAFFOLDS);
				}else if(a.startsWith("maxlen") || a.startsWith("chromlen")){
					long len=Parse.parseKMG(b);
					assert(len>0 && len<=Integer.MAX_VALUE);
					MAX_LENGTH=(int)len;
				}else if(a.equals("writechroms")){
					writeChroms=Parse.parseBoolean(b);
				}else if(a.equals("chromgz") || a.equals("gz")){
					Data.CHROMGZ=Parse.parseBoolean(b);
				}else if(a.equals("retain")){
					RETAIN=Parse.parseBoolean(b);
				}else if(a.equals("scafprefixes")){
					scafprefixes=Parse.parseBoolean(b);
				}else if(a.equals("ziplevel") || a.equals("zl")){
					ReadWrite.ZIPLEVEL=Integer.parseInt(b);
				}else if(a.equals("waitforwriting")){
					WAIT_FOR_WRITING=Parse.parseBoolean(b);
				}else{
					if(i>2){
						System.err.println("Unknown parameter "+args[i]);
//						throw new RuntimeException("Unknown parameter "+args[i]);
					}
				}
			}
		}
		
		WAIT_FOR_WRITING=(WAIT_FOR_WRITING || ReadWrite.USE_GZIP || ReadWrite.USE_PIGZ);
		
		ArrayList<ChromosomeArray> r=RETAIN ? new ArrayList<ChromosomeArray>() : null;
		
//		assert(false) : Arrays.toString(args);
//		assert(RETAIN);
		
		if(genome<0){genome=Integer.parseInt(args[1]);} //Legacy
		if(genome<0){throw new RuntimeException("Please specify a genome build number.");}
		
		if(writeinfo){
			if(chroms<0){chroms=Integer.parseInt(args[2]);} //Legacy
			if(chroms<0){throw new RuntimeException("Please the number of chroms.");}
			writeInfo(genome, chroms, name, null, false, scafprefixes);
		}else{
			if(infile==null){infile=args[0].replace('\\', '/');} //Legacy
			if(infile==null){throw new RuntimeException("Please specify an input file.");}
			{
				File f=new File(infile);
				if(!f.exists() || f.isDirectory()){
					if(!infile.startsWith("stdin")){
						throw new RuntimeException("Not a valid file: "+f);
					}
				}
			}
			String outRoot=Data.ROOT_GENOME+genome+"/";
			
			FastaToChromArrays2 ftca=new FastaToChromArrays2();
			ftca.makeChroms(infile, outRoot, name, genScaffoldInfo, writeChroms, r, scafprefixes);
		}
		
		WRITE_IN_THREAD=oldWIT;
		return r;
	}
	
	private FastaToChromArrays2(){}
	
	
	private static int[] countInfo(ChromosomeArray ca){
		int contigs=0;
		int startPad=0;
		int stopPad=0;
		int undefined=0;
		int defined=0;//=ca.countDefinedBases();
		
		int lastN=-1;
		int lastDef=-1;
		
		for(int i=0; i<=ca.maxIndex; i++){
			byte b=ca.get(i);
			if(AminoAcid.isFullyDefined(b)){
				if(defined==0){startPad=i; contigs++;}
				else if(i-lastDef>contigTrigger){contigs++;}
				lastDef=i;
				defined++;
			}else{
				lastN=i;
				undefined++;
			}
		}
		
		if(contigs>0 && lastN==ca.maxIndex){
			stopPad=lastN-lastDef;
		}else{
//			System.err.println(lastN+", "+lastDef+", "+ca.maxIndex);
		}
		
		return new int[] {ca.chromosome, 1, contigs, (ca.maxIndex+1), defined, undefined, startPad, stopPad};
	}
	
	@Deprecated
	public static void writeInfo(int genome, int chroms, String name, String source, boolean unload, boolean scafNamePrefix){
		Data.GENOME_BUILD=genome;
		Data.chromosomePlusMatrix=new ChromosomeArray[chroms+1];
				
		String outRoot=Data.ROOT_GENOME+genome+"/";
		TextStreamWriter info=new TextStreamWriter(outRoot+"info.txt", true, false, false);
		info.start();
		info.print("#Chromosome sizes\n");
		try {
			info.print("#Generated on\t"+new Date()+"\n");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		info.print("#Version\t"+VERSION+"\n");
		info.print("#chrom\tscaffolds\tcontigs\tlength\tdefined\tundefined\tstartPad\tstopPad\n");
		
		
		long bases=0;
		long definedBases=0;
		
		long contigSum=0;
		
		for(int chrom=1; chrom<=chroms; chrom++){
			ChromosomeArray ca=Data.getChromosome(chrom);
			int[] v=countInfo(ca);
			info.print(v[0]+"\t"+v[1]+"\t"+v[2]+"\t"+v[3]+"\t"+v[4]+"\t"+v[5]+"\t"+v[6]+"\t"+v[7]+"\n");
			
			bases+=v[3];
			definedBases+=v[4];
			contigSum+=v[2];
			if(unload){Data.unload(chrom, false);}
		}
		info.poison();
		StringBuilder sb=new StringBuilder();
		sb.append("#Summary\n");
		try {
			sb.append("#Generated on\t"+new Date()+"\n");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		sb.append("#Version\t"+VERSION+"\n");
		sb.append("chroms\t"+(chroms)+"\n");
		sb.append("bases\t"+bases+"\n");
		sb.append("defined\t"+definedBases+"\n");
		sb.append("undefined\t"+(bases-definedBases)+"\n");
		sb.append("contigs\t"+contigSum+"\n");
		sb.append("scaffolds\t"+chroms+"\n");
		sb.append("interpad\t"+MID_PADDING+"\n");
		if(name!=null){sb.append("name\t"+name+"\n");}
		if(source!=null){sb.append("source\t"+source+"\n");}
		if(scafNamePrefix){sb.append("scafprefixes\t"+scafNamePrefix+"\n");}//else{assert(false);}
		ReadWrite.writeString(sb, outRoot+"summary.txt", false);
		info.waitForFinish();
	}
	
	private int makeChroms(String fname, String outRoot, String genomeName, boolean genScaffolds, boolean writeChroms, ArrayList<ChromosomeArray> r,
			boolean scafNamePrefix){
		
		if(!NODISK){
			File f=new File(outRoot);
			if(!f.exists()){
				if(!NODISK){f.mkdirs();}
			}else if(overwrite){
				for(File g : f.listFiles()){
					String s=g.getName();
					if(g.isFile() && s.contains(".chrom")){
						System.err.println("Deleting "+s);
						g.delete();
					}
				}
			}
			
			f=new File(outRoot.replace("ref/genome/", "ref/index/"));
			if(!f.exists()){
				if(!NODISK){f.mkdirs();}
			}else if(overwrite){
				for(File g : f.listFiles()){
					String s=g.getName();
					if(g.isFile() && (s.endsWith(".int2d") || s.endsWith(".block") || s.endsWith(".block2.gz") || s.endsWith(".blockB") || s.endsWith(".blockB2.gz"))){
						System.err.println("Deleting "+s);
						g.delete();
					}
				}
			}
		}
		
		ByteFile1 tf=new ByteFile1(fname, false);
		int chrom=1;
		
		TextStreamWriter infoWriter=null, scafWriter=null;
		ArrayList<String> infolist=null, scaflist=null;

		if(NODISK){
			infolist=new ArrayList<String>();
			infolist.add("#Chromosome sizes");
			try {
				infolist.add("#Generated on\t"+new Date());
			} catch (Exception e1) {
				e1.printStackTrace();
			}
			infolist.add("#Version\t"+VERSION);
			infolist.add("#chrom\tscaffolds\tcontigs\tlength\tdefined\tundefined\tstartPad\tstopPad");
		}else{
			infoWriter=new TextStreamWriter(outRoot+"info.txt", true, false, false);
			infoWriter.start();
			infoWriter.print("#Chromosome sizes\n");
			try {
				//			System.err.println(new Date());
				infoWriter.print("#Generated on\t"+new Date()+"\n");
			} catch (Exception e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			infoWriter.print("#Version\t"+VERSION+"\n");
			infoWriter.print("#chrom\tscaffolds\tcontigs\tlength\tdefined\tundefined\tstartPad\tstopPad\n");
		}
		
		if(genScaffolds){
			if(NODISK){
				scaflist=new ArrayList<String>();
				scaflist.add("#Scaffold names");
				try {
					scaflist.add("#Generated on\t"+new Date());
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				scaflist.add("#Version\t"+VERSION);
				scaflist.add("#chrom\tid\tstart\tlength\tname");
			}else{
				//System.err.println("*123 Making ScafWriter; "+ReadWrite.countActiveThreads()+", "+ReadWrite.USE_GZIP+", "+ReadWrite.USE_PIGZ);
				scafWriter=new TextStreamWriter(outRoot+"scaffolds.txt.gz", true, false, false);
				scafWriter.start();
				scafWriter.print("#Scaffold names\n");
				try {
					scafWriter.print("#Generated on\t"+new Date()+"\n");
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				scafWriter.print("#Version\t"+VERSION+"\n");
				scafWriter.print("#chrom\tid\tstart\tlength\tname\n");
			}
		}
		
		
		for(ChromosomeArray ca=makeNextChrom(tf, chrom, infoWriter, scafWriter, infolist, scaflist); ca!=null;
				ca=makeNextChrom(tf, chrom, infoWriter, scafWriter, infolist, scaflist)){
			if(ca.array.length>ca.maxIndex+1){ca.resize(ca.maxIndex+1);}
			if(RETAIN){r.add(ca);}
			
			if(writeChroms){
				String x=outRoot+"chr"+chrom+Data.chromExtension();
				if(new File(x).exists() && !overwrite){throw new RuntimeException("Tried to overwrite existing file "+x+", but overwrite=false.");}
				ReadWrite.writeObjectInThread(ca, x, false);
				System.err.println("Writing chunk "+chrom);
			}
			chrom++;
		}
		lastHeader=nextHeader=null;
		
		tf.close();
		if(infoWriter!=null){infoWriter.poison();}
		if(scafWriter!=null){
			//System.err.println("*123 Killing ScafWriter; "+ReadWrite.countActiveThreads()+", "+ReadWrite.USE_GZIP+", "+ReadWrite.USE_PIGZ);
			scafWriter.poison();
		}
		
		StringBuilder sb=new StringBuilder();
		sb.append("#Summary\n");
		try {
			sb.append("#Generated on\t"+new Date()+"\n");
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		sb.append("#Version\t"+VERSION+"\n");
		sb.append("chroms\t"+(chrom-1)+"\n");
		sb.append("bases\t"+lengthSum+"\n");
		assert((definedSum+undefinedSum)==lengthSum) : definedSum+", "+undefinedSum+", "+lengthSum;
		sb.append("defined\t"+definedSum+"\n");
		sb.append("undefined\t"+undefinedSum+"\n");
		sb.append("contigs\t"+contigSum+"\n");
		sb.append("scaffolds\t"+scaffoldSum+"\n");
		sb.append("interpad\t"+MID_PADDING+"\n");
		if(genomeName!=null){sb.append("name\t"+genomeName+"\n");}
		if(fname!=null){
			File f=new File(fname);
			String cpath=null;
			try {
				cpath=f.getCanonicalPath();
			} catch (IOException e) {
				cpath=f.getAbsolutePath();
			}
			sb.append("source\t"+cpath+"\n");
			sb.append("bytes\t"+f.length()+"\n");
			sb.append("last modified\t"+f.lastModified()+"\n");
		}
		if(scafNamePrefix){sb.append("scafprefixes\t"+scafNamePrefix+"\n");}//else{assert(false);}
		if(NODISK){
			SUMMARY_LIST=new ArrayList<String>();
			String[] split=sb.toString().split("\n");
			for(String s : split){SUMMARY_LIST.add(s);}
		}else{
			ReadWrite.writeString(sb, outRoot+"summary.txt", false);
		}
		
		if(infoWriter!=null){infoWriter.waitForFinish();}
		if(scafWriter!=null){
			//System.err.println("*123 Waiting For ScafWriter; "+ReadWrite.countActiveThreads()+", "+ReadWrite.USE_GZIP+", "+ReadWrite.USE_PIGZ);
			scafWriter.waitForFinish();
			//System.err.println("*123 ScafWriter Finished; "+ReadWrite.countActiveThreads()+", "+ReadWrite.USE_GZIP+", "+ReadWrite.USE_PIGZ);
		}
		
		if(WAIT_FOR_WRITING && ReadWrite.countActiveThreads()>0){
			System.err.println("Waiting for writing to finish.");
			ReadWrite.waitForWritingToFinish();
			System.err.println("Finished.");
			//System.err.println("*123 countActiveThreads Finished; "+ReadWrite.countActiveThreads()+", "+ReadWrite.USE_GZIP+", "+ReadWrite.USE_PIGZ);
		}

		if(infolist!=null){
			INFO_LIST=infolist;
			LISTBUILD=Data.GENOME_BUILD;
		}else{INFO_LIST=null;}
		if(scaflist!=null){
			SCAF_LIST=scaflist;
			LISTBUILD=Data.GENOME_BUILD;
		}else{SCAF_LIST=null;}
		
		return chrom-1;
	}
	
	private ChromosomeArray makeNextChrom(ByteFile1 tf, int chrom, TextStreamWriter infoWriter, TextStreamWriter scafWriter, ArrayList<String> infolist, ArrayList<String> scaflist){
		ChromosomeArray ca=new ChromosomeArray(chrom, (byte)Shared.PLUS, 0, 120000+START_PADDING);
		ca.maxIndex=-1;
		for(int i=0; i<START_PADDING; i++){ca.set(i, 'N');}
		
		if(verbose){System.err.println("chrom="+chrom+", lastHeader="+lastHeader+", nextHeader="+nextHeader);}
		
		int scaffolds=0;
		if(currentScaffold!=null && currentScaffold.length()>0){
			assert(currentScaffold.length()>0) : currentScaffold.length();
			assert(lastHeader!=null);
			assert(currentScaffold.length()+END_PADDING+ca.maxIndex<MAX_LENGTH) : currentScaffold.length()+", "+END_PADDING+", "+ca.maxIndex+", "+MAX_LENGTH;
			
//			System.err.println("A: Writing a scaffold because currentScaffold = "+currentScaffold);
			scaffoldSum++;
			if(scafWriter!=null){scafWriter.print(chrom+"\t"+scaffoldSum+"\t"+(ca.maxIndex+1)+"\t"+currentScaffold.length()+"\t"+lastHeader+"\n");}
			if(scaflist!=null && lastHeader!=null){
				scaflist.add(chrom+"\t"+scaffoldSum+"\t"+(ca.maxIndex+1)+"\t"+currentScaffold.length()+"\t"+lastHeader);
				if(verbose){System.err.println("A: Added to scaflist: "+scaflist.get(scaflist.size()-1));}
			}
			ca.set(ca.maxIndex+1, currentScaffold);
			scaffolds++;
			
			currentScaffold.setLength(0);
			lastHeader=nextHeader;
		}
		
//		if()
		
		while((currentScaffold=nextScaffold(currentScaffold, tf))!=null){
			if(currentScaffold.length()+MID_PADDING+END_PADDING+ca.maxIndex>MAX_LENGTH){break;}
			if(scaffolds>0 && !MERGE_SCAFFOLDS){break;}
			
			if(scaffolds>0){
				for(int i=0; i<MID_PADDING; i++){
					ca.set(ca.maxIndex+1, 'N');
				}
			}
			if(currentScaffold.length()>=MIN_SCAFFOLD){
//				System.err.println("B: Writing a scaffold because currentScaffold = "+currentScaffold);
				scaffoldSum++;
				if(scafWriter!=null){scafWriter.print(chrom+"\t"+scaffoldSum+"\t"+(ca.maxIndex+1)+"\t"+currentScaffold.length()+"\t"+lastHeader+"\n");}
				if(scaflist!=null){
					scaflist.add(chrom+"\t"+scaffoldSum+"\t"+(ca.maxIndex+1)+"\t"+currentScaffold.length()+"\t"+lastHeader);
					if(verbose){System.err.println("B: Added to scaflist: "+scaflist.get(scaflist.size()-1));}
				}
				ca.set(ca.maxIndex+1, currentScaffold);
				scaffolds++;
			}
			
			currentScaffold.setLength(0);
			lastHeader=nextHeader;
		}
		
		if(verbose){System.err.println("lastHeader="+lastHeader);}
		
		if(scaffolds==0){return null;}
		
		if(END_PADDING>0){
			int terminalN=0;
			for(int i=ca.maxIndex; i>=0 && terminalN<END_PADDING; i--){
				if(ca.get(i)=='N'){terminalN++;}
				else{break;}
			}
//			System.err.println("\nAdding Ns: ref.length="+ca.maxIndex);
			while(terminalN<=END_PADDING && ca.maxIndex<MAX_LENGTH-1){
//				System.out.print("N");
				ca.set(ca.maxIndex+1, 'N');
				terminalN++;
			}
//			System.err.println("\nAdded Ns: ref.length="+ca.maxIndex);
		}
		
		int[] v=countInfo(ca);
		v[6]=Tools.max(0, Tools.min(START_PADDING, v[6])); //In case input scaffolds had leading undefined bases
		v[7]=Tools.max(0, Tools.min(END_PADDING, v[7])); //In case input scaffolds had trailing undefined bases
		if(infoWriter!=null){
//			infoWriter.print("#chrom\tscaffolds\tcontigs\tlength\tdefined\tundefined\tstartPad\tstopPad\n");
			infoWriter.print(v[0]+"\t"+scaffolds+"\t"+v[2]+"\t"+v[3]+"\t"+v[4]+"\t"+v[5]+"\t"+v[6]+"\t"+v[7]+"\n");
		}
		if(infolist!=null){
			infolist.add(v[0]+"\t"+scaffolds+"\t"+v[2]+"\t"+v[3]+"\t"+v[4]+"\t"+v[5]+"\t"+v[6]+"\t"+v[7]);
		}
		lengthSum+=v[3];
		definedSum+=v[4];
		undefinedSum+=v[5];
		contigSum+=v[2];

		assert((definedSum+undefinedSum)==lengthSum) : definedSum+", "+undefinedSum+", "+lengthSum+
			"; "+ca.countDefinedBases()+", "+(ca.maxIndex+1)+"\n"+ca.getString(0, ca.maxIndex);
		
		return ca;
	}
	
	
	private ByteBuilder nextScaffold(ByteBuilder sb, ByteFile1 tf){
		if(sb==null){sb=new ByteBuilder(100);}
		else{sb.setLength(0);}
		if(!tf.isOpen()){return null;}
		
		byte[] s=tf.nextLine();
		
		for(; s!=null && (s.length==0 || s[0]!='>'); s=tf.nextLine()){
//			if(TRANSLATE_U_TO_T){
//				for(int i=0; i<s.length; i++){
//					if(s[i]=='U'){s[i]='T';}
//				}
//			}
			sb.append(s);
//			for(int i=0; i<s.length(); i++){
//				char c=s.charAt(i);
//				assert(Tools.isLetter(c));
//				sb.append(Tools.toUpperCase(c));
//			}
		}
		
		nextHeader=(s==null ? null : new String(s, 1, s.length-1));
		if(s==null && sb.length()==0){return null;}
		return sb;
	}
	
	private String lastHeader;
	private String nextHeader;
	private ByteBuilder currentScaffold;
	private long scaffoldSum=0;
	private long lengthSum=0;
	private long definedSum=0;
	private long undefinedSum=0;
	private long contigSum=0;
	
	
	public static final int currentVersion(){return VERSION;}
	
	public static boolean MERGE_SCAFFOLDS=true;
	public static boolean WRITE_IN_THREAD=false;
	public static boolean overwrite=true;
	public static boolean append=false;
	public static int START_PADDING=8000; //Always applied
	public static int MID_PADDING=300; //Applied when merging scaffolds
	public static int END_PADDING=8000; //Only applied if not enough terminal Ns
	public static int MIN_SCAFFOLD=1;
	public static int contigTrigger=10;
	public static int VERSION=5;
	public static int MAX_LENGTH=(1<<29)-200000;
//	public static boolean TRANSLATE_U_TO_T;
	
	public static boolean verbose=false;
	public static boolean RETAIN=false;
	public static boolean WAIT_FOR_WRITING=true;
	public static boolean NODISK=false;
	public static int LISTBUILD=-1;
	public static ArrayList<String> INFO_LIST, SCAF_LIST, SUMMARY_LIST;
	
//	public static boolean GENERATE_SCAFFOLD_INFO=true;
	
}
