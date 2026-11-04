package dna;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.CrisWrapper;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;

/**
 * Legacy FASTA/FASTQ-to-chromosome converter using a read-stream wrapper.
 * Copies sequence bases into chromosome arrays; input qualities are not retained.
 * Configuration and published metadata are shared static state, so callers must
 * serialize use. Several temporary input settings are restored only on normal return.
 *
 * Note: this class currently has NO callers tree-wide (no .sh, no importers); the live
 * canonical reference indexer is {@link FastaToChromArrays2} (used by BBMap5, RefToIndex,
 * Data, AssemblyStats2, etc.). Its makeNextChrom is the un-hardened predecessor of
 * FastaToChromArrays2's and lacks that class's later oversize-scaffold guards (its #001
 * KillSwitch.kill on scaffolds &gt; MAX_SINGLE_SCAFFOLD, #002 minscaf on carried scaffolds,
 * #003 single-large-scaffold-own-chrom): here the oversize check is a bare assert.
 * An assertion can leave already-started writers waiting; there is no KillSwitch
 * or MAX_SINGLE_SCAFFOLD handling here. Preserve this legacy implementation separately
 * from the live converter rather than assuming equivalent boundary behavior.
 * @author Brian Bushnell
 * @date Jul 18, 2014
 *
 */
public class ChromArrayMaker{

//	Example:
//	dna.ChromArrayMaker ecoli_K12.fa 1 writeinthread=false genscaffoldinfo=true retain waitforwriting=false
//	gzip=true chromc=false maxlen=536670912 writechroms=true minscaf=1 midpad=300 startpad=8000 stoppad=8000 nodisk=false

	/** Program entry point that delegates to main2.
	 * @param args Command-line arguments for chromosome array creation */
	public static void main(String[] args){
		main2(args);
	}

	/**
	 * Main processing method that parses arguments and creates chromosome arrays.
	 * Supports named options and legacy positional input/build/chromosome-count arguments.
	 * NODISK redirects metadata into static lists; writechroms independently controls
	 * chromosome-object output. RETAIN controls the returned list, including an empty
	 * list in writeinfo mode. FASTQ pairing detection is disabled during conversion.
	 * Shared path/compression/conversion options persist after the call; WRITE_IN_THREAD
	 * and pairing flags are restored only when control reaches their normal cleanup.
	 *
	 * @param args Command-line arguments including input file, genome build, and options
	 * @return ArrayList of ChromosomeArray objects if RETAIN=true, null otherwise
	 */
	public static ArrayList<ChromosomeArray> main2(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null, false);
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

			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(Parser.parseZip(arg, a, b)){
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
			}else if(a.equals("waitforwriting")){
				WAIT_FOR_WRITING=Parse.parseBoolean(b);
			}else{
				if(i>2){
					System.err.println("WARNING: Unknown parameter "+args[i]);
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

			ChromArrayMaker ftca=new ChromArrayMaker();

			{
				//TODO: Probable bug (DNA-018) - exceptions skip restoration of pairing and
				//other temporary input settings. No production callers found; repair needs a bounded failure test.
				boolean oldTI=FASTQ.TEST_INTERLEAVED;
				boolean oldFI=FASTQ.FORCE_INTERLEAVED;
				FASTQ.TEST_INTERLEAVED=false;
				FASTQ.FORCE_INTERLEAVED=false;
				ftca.makeChroms(infile, outRoot, name, genScaffoldInfo, writeChroms, r, scafprefixes);
				FASTQ.TEST_INTERLEAVED=oldTI;
				FASTQ.FORCE_INTERLEAVED=oldFI;
			}
		}

		WRITE_IN_THREAD=oldWIT;
		return r;
	}

	/** Private constructor to prevent external instantiation */
	private ChromArrayMaker(){}

	/**
	 * Analyzes a chromosome array to count contigs, bases, and padding information.
	 * Scans coordinates zero through the inclusive maximum. A new contig starts when
	 * consecutive defined-base positions differ by more than contigTrigger. Fully
	 * undefined arrays have zero contigs and zero reported terminal padding.
	 *
	 * @param ca The chromosome array to analyze
	 * @return Array containing [chromosome, 1, contigs, length, defined,
	 * undefined, startPad, stopPad]
	 */
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

	/**
	 * Writes chromosome information and summary files for an existing genome build.
	 * Replaces Data's genome build and chromosome cache, then reads existing arrays.
	 * Each chromosome is counted as one scaffold. Joins the info writer on normal
	 * completion, but currently does not check its retained error state.
	 *
	 * @param genome Genome build number
	 * @param chroms Number of chromosomes to process
	 * @param name Organism name for summary file
	 * @param source Source file path for summary file
	 * @param unload Whether to unload chromosome arrays after processing
	 * @param scafNamePrefix Value recorded as the summary's scaffold-prefix flag
	 * @deprecated Use makeChroms method instead for new conversions
	 */
	@Deprecated
	public static void writeInfo(int genome, int chroms, String name, String source, boolean unload, boolean scafNamePrefix){
		Data.GENOME_BUILD=genome;
		Data.chromosomePlusMatrix=new ChromosomeArray[chroms+1];

		String outRoot=Data.ROOT_GENOME+genome+"/";
		TextStreamWriter info=new TextStreamWriter(outRoot+"info.txt", true, false, false);
		info.start();
		info.print("#Chromosome sizes\n");
		try{
			info.print("#Generated on\t"+new Date()+"\n");
		}catch(Exception e){
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
		try{
			sb.append("#Generated on\t"+new Date()+"\n");
		}catch(Exception e){
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
		//TODO: Probable bug (DNA-017) - joining does not propagate info.errorState.
	}

	/**
	 * Creates chromosome arrays from input FASTA/FASTQ file with configurable options.
	 * Reads sequences, applies padding, merges scaffolds based on settings, and writes output files.
	 * Metadata goes to files or static lists according to NODISK; writeChroms separately
	 * requests background serialization. With overwrite enabled, existing chromosome
	 * and selected index files are removed before conversion. The input wrapper closes
	 * itself at exhaustion. Metadata writers are joined on normal completion, but
	 * retained reader/writer errors are not currently checked before publishing summary.
	 *
	 * @param fname Input FASTA or FASTQ file path
	 * @param outRoot Output directory path ending in a separator
	 * @param genomeName Organism name for metadata files
	 * @param genScaffolds Whether to generate scaffold information file
	 * @param writeChroms Whether to write individual chromosome array files
	 * @param r Mutable destination list, required when RETAIN=true
	 * @param scafNamePrefix Flag recorded in summary metadata; does not rewrite read IDs
	 * @return Number of chromosomes created
	 */
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

		final FileFormat ffin=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);

		final boolean OLD_SPLIT_READS=FastaReadInputStream.SPLIT_READS;
		FastaReadInputStream.SPLIT_READS=false;
		final int oldNum=Shared.numBuffers();
		Shared.setBuffers(4);
		final CrisWrapper criswrapper=new CrisWrapper(-1, false, ffin, null, null, null);
		Shared.setBuffers(oldNum);

		int chrom=1;

		TextStreamWriter infoWriter=null, scafWriter=null;
		ArrayList<String> infolist=null, scaflist=null;

		if(NODISK){
			infolist=new ArrayList<String>();
			infolist.add("#Chromosome sizes");
			try{
				infolist.add("#Generated on\t"+new Date());
			}catch(Exception e1){
				e1.printStackTrace();
			}
			infolist.add("#Version\t"+VERSION);
			infolist.add("#chrom\tscaffolds\tcontigs\tlength\tdefined\tundefined\tstartPad\tstopPad");
		}else{
			infoWriter=new TextStreamWriter(outRoot+"info.txt", true, false, false);
			infoWriter.start();
			infoWriter.print("#Chromosome sizes\n");
			try{
				//			System.err.println(new Date());
				infoWriter.print("#Generated on\t"+new Date()+"\n");
			}catch(Exception e1){
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
				try{
					scaflist.add("#Generated on\t"+new Date());
				}catch(Exception e){
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
				try{
					scafWriter.print("#Generated on\t"+new Date()+"\n");
				}catch(Exception e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				scafWriter.print("#Version\t"+VERSION+"\n");
				scafWriter.print("#chrom\tid\tstart\tlength\tname\n");
			}
		}

		for(ChromosomeArray ca=makeNextChrom(criswrapper, chrom, infoWriter, scafWriter, infolist, scaflist); ca!=null;
				ca=makeNextChrom(criswrapper, chrom, infoWriter, scafWriter, infolist, scaflist)){
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

		//TODO: Probable bug (DNA-017) - CrisWrapper.errorState and the metadata writers'
		//errorState are ignored, so close failures can still publish a completion summary.
		//No production caller found; the live FastaToChromArrays2 repair is separate.
		FastaReadInputStream.SPLIT_READS=OLD_SPLIT_READS;

		if(infoWriter!=null){infoWriter.poison();}
		if(scafWriter!=null){
			//System.err.println("*123 Killing ScafWriter; "+ReadWrite.countActiveThreads()+", "+ReadWrite.USE_GZIP+", "+ReadWrite.USE_PIGZ);
			scafWriter.poison();
		}

		StringBuilder sb=new StringBuilder();
		sb.append("#Summary\n");
		try{
			sb.append("#Generated on\t"+new Date()+"\n");
		}catch(Exception e1){
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
			try{
				cpath=f.getCanonicalPath();
			}catch(IOException e){
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

	/**
	 * Creates a single chromosome array by reading and merging scaffolds from input.
	 * Applies start padding, merges scaffolds based on length limits and settings,
	 * and applies end padding.
	 * Updates cumulative statistics and writes scaffold information. currentScaffold
	 * carries one unconsumed read across calls when a packing boundary is reached.
	 * The carried-read path does not apply MIN_SCAFFOLD, and a first oversized read
	 * can return null before it is consumed; these are legacy limitations, not guarantees.
	 *
	 * @param criswrapper Input stream wrapper for reading sequences
	 * @param chrom Chromosome number being created
	 * @param infoWriter Writer for info.txt file (null if NODISK=true)
	 * @param scafWriter Writer for scaffolds.txt.gz file (null if not generating scaffold info)
	 * @param infolist List to store info lines if NODISK=true
	 * @param scaflist List to store scaffold lines if NODISK=true
	 * @return ChromosomeArray containing merged scaffolds, or null if no scaffolds processed
	 */
	private ChromosomeArray makeNextChrom(CrisWrapper criswrapper, int chrom, TextStreamWriter infoWriter, TextStreamWriter scafWriter, ArrayList<String> infolist, ArrayList<String> scaflist){
		assert(FastaReadInputStream.SPLIT_READS==false);
		ChromosomeArray ca=new ChromosomeArray(chrom, (byte)Shared.PLUS, 0, 120000+START_PADDING);
		ca.maxIndex=-1;
		for(int i=0; i<START_PADDING; i++){ca.set(i, 'N');}

		if(verbose){System.err.println("chrom="+chrom+", currentScaffold="+(currentScaffold==null ? "null" : currentScaffold.id));}

		int scaffolds=0;
		if(currentScaffold!=null && currentScaffold.length()>0){
			assert(currentScaffold.length()>0);
			assert(currentScaffold.length()+END_PADDING+ca.maxIndex<MAX_LENGTH);

			if(verbose){System.err.println("A: Writing a scaffold because currentScaffold = "+currentScaffold);}
			scaffoldSum++;
			if(scafWriter!=null){scafWriter.print(chrom+"\t"+scaffoldSum+"\t"+(ca.maxIndex+1)+"\t"+currentScaffold.length()+"\t"+currentScaffold.id+"\n");}
			if(scaflist!=null && currentScaffold!=null){
				scaflist.add(chrom+"\t"+scaffoldSum+"\t"+(ca.maxIndex+1)+"\t"+currentScaffold.length()+"\t"+currentScaffold.id);
				if(verbose){System.err.println("A: Added to scaflist: "+scaflist.get(scaflist.size()-1));}
			}
			ca.set(ca.maxIndex+1, currentScaffold.bases);
			scaffolds++;

			currentScaffold=null;
		}

		while((currentScaffold=criswrapper.next())!=null){
			if(currentScaffold.length()+MID_PADDING+END_PADDING+ca.maxIndex>MAX_LENGTH){break;}
			if(scaffolds>0 && !MERGE_SCAFFOLDS){break;}

			if(scaffolds>0){
				for(int i=0; i<MID_PADDING; i++){
					ca.set(ca.maxIndex+1, 'N');
				}
			}
			if(currentScaffold.length()>=MIN_SCAFFOLD){
				if(verbose){System.err.println("B: Writing a scaffold because currentScaffold = "+currentScaffold);}
				scaffoldSum++;
				if(scafWriter!=null){scafWriter.print(chrom+"\t"+scaffoldSum+"\t"+(ca.maxIndex+1)+"\t"+currentScaffold.length()+"\t"+currentScaffold.id+"\n");}
				if(scaflist!=null){
					scaflist.add(chrom+"\t"+scaffoldSum+"\t"+(ca.maxIndex+1)+"\t"+currentScaffold.length()+"\t"+currentScaffold.id);
					if(verbose){System.err.println("B: Added to scaflist: "+scaflist.get(scaflist.size()-1));}
				}
				ca.set(ca.maxIndex+1, currentScaffold.bases);
				scaffolds++;
			}
			currentScaffold=null;
		}

//		if(verbose){System.err.println("lastHeader="+lastHeader);}

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

//	private String lastHeader;
//	private String nextHeader;
	/** Pending input read, retained across chromosome boundaries; only its bases and ID are used. */
	private Read currentScaffold;
	/** Running total of scaffolds processed across all chromosomes */
	private long scaffoldSum=0;
	/** Running total of bases across all chromosomes */
	private long lengthSum=0;
	/** Running total of defined (non-N) bases across all chromosomes */
	private long definedSum=0;
	/** Running total of undefined (N) bases across all chromosomes */
	private long undefinedSum=0;
	/** Running total of contigs across all chromosomes */
	private long contigSum=0;

	/** Returns the current version number for chromosome array format */
	public static final int currentVersion(){return VERSION;}

	/** Whether to merge multiple scaffolds into single chromosomes */
	public static boolean MERGE_SCAFFOLDS=true;
	/** Legacy parsed option; chromosome output uses the background writer regardless of this value. */
	public static boolean WRITE_IN_THREAD=false;
	/** Whether to overwrite existing output files */
	public static boolean overwrite=true;
	/** Legacy parsed option; not consulted by this implementation's output calls. */
	public static boolean append=false;
	/** Number of N bases to add at the start of each chromosome */
	public static int START_PADDING=8000; //Always applied
	/** Number of N bases to insert between merged scaffolds */
	public static int MID_PADDING=300; //Applied when merging scaffolds
	/** Legacy terminal-padding threshold; extension uses an inclusive comparison and a length cap. */
	public static int END_PADDING=8000;
	/** Minimum length for newly read scaffolds; the carried-over path currently bypasses this filter. */
	public static int MIN_SCAFFOLD=1;
	/** Gap size threshold for counting separate contigs within scaffolds */
	public static int contigTrigger=10;
	/** Version number for chromosome array format compatibility */
	public static int VERSION=5;
	/** Legacy packing threshold; not an unconditional bound for every input/parameter combination. */
	public static int MAX_LENGTH=(1<<29)-200000;

	/** Whether to print detailed processing information */
	public static boolean verbose=false;
	/** Whether to keep chromosome arrays in memory after creation */
	public static boolean RETAIN=false;
	/**
	 * Whether to wait for the shared background object-writer counter after metadata output.
	 * main2 also enables this when gzip or pigz output is enabled.
	 */
	public static boolean WAIT_FOR_WRITING=true;
	/** Redirects metadata to memory and skips directory cleanup; writeChroms remains independent. */
	public static boolean NODISK=false;
	/** Genome build number associated with the in-memory lists */
	public static int LISTBUILD=-1;
	/** Shared metadata destinations; updated on normal completion and not thread-safe. */
	public static ArrayList<String> INFO_LIST, SCAF_LIST, SUMMARY_LIST;

//	public static boolean GENERATE_SCAFFOLD_INFO=true;

}
