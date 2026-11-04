package dna;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

import fileIO.ByteFile1;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Builds chromosome arrays and reference metadata from byte-oriented FASTA input.
 * This is the reference-conversion path used by RefToIndex. It packs scaffolds
 * using configurable padding/length limits and can retain arrays for its caller.
 * Configuration, Data paths and published metadata lists are shared static state;
 * callers must serialize conversions and arrange the desired settings each time.
 *
 * Metadata completion checks reject retained input or writer errors before a new
 * summary is published. This is not a transactional build: earlier chromosome or
 * metadata files may already exist after failure, and prior static lists are not
 * cleared at entry. The deprecated writeInfo path has separate legacy error handling.
 * @author Brian Bushnell
 * @date Jul 30, 2014
 *
 */
public class FastaToChromArrays2{

//	Example:
//	jgi.FastaToChromArrays ecoli_K12.fa 1 writeinthread=false genscaffoldinfo=true retain waitforwriting=false
//	gzip=true chromc=false maxlen=536670912 writechroms=true minscaf=1 midpad=300 startpad=8000 stoppad=8000 nodisk=false

	/** Program entry point that delegates to main2.
	 * @param args Command-line arguments for FASTA processing */
	public static void main(String[] args){
		main2(args);
	}

	/**
	 * Main processing method that parses arguments and converts FASTA to ChromosomeArrays.
	 * Handles both info-only mode and full conversion mode.
	 * Supports named options and legacy positional input/build/chromosome-count values.
	 * NODISK redirects metadata to static lists; writechroms independently controls
	 * chromosome serialization. RETAIN selects the returned list, which is empty in
	 * writeinfo mode. Other shared options persist; the retained legacy WRITE_IN_THREAD
	 * state is restored on normal and exceptional exits after preprocessing. This does
	 * not control chromosome-write scheduling or provide concurrent/transactional reset.
	 *
	 * @param args Command-line arguments including input file, genome build, and options
	 * @return ArrayList of ChromosomeArray objects if RETAIN=true, null otherwise
	 * @throws RuntimeException If input validation or checked reference completion fails
	 */
	public static ArrayList<ChromosomeArray> main2(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}

		boolean oldWIT=WRITE_IN_THREAD;
		try{
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

			return r;
		}finally{
			WRITE_IN_THREAD=oldWIT;//DNA-023: restore only the temporary legacy state.
		}
	}

	/** Private constructor to prevent external instantiation */
	private FastaToChromArrays2(){}

	/**
	 * Counts statistics for a ChromosomeArray including contigs, defined/undefined bases, and padding.
	 * Identifies contig boundaries when undefined gaps contain at least contigTrigger bases.
	 * Scans coordinates zero through the inclusive maximum. Contigs split when the
	 * distance between defined positions exceeds contigTrigger. Fully undefined
	 * arrays report zero contigs and zero terminal padding.
	 *
	 * @param ca The ChromosomeArray to analyze
	 * @return Array: [chromosome, scaffolds=1, contigs, length, defined, undefined, startPad, stopPad]
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

		return new int[]{ca.chromosome, 1, contigs, (ca.maxIndex+1), defined, undefined, startPad, stopPad};
	}

	/**
	 * Writes genome information files for pre-existing ChromosomeArrays.
	 * Creates info.txt with per-chromosome statistics and summary.txt with totals.
	 * Replaces Data's build and chromosome cache and treats each chromosome as one
	 * scaffold. Joins a successfully started info writer on normal and exceptional
	 * producer exits. Preserves the original unchecked failure, suppressing distinct
	 * cleanup failures; retained writer errors reject before writing a new summary.
	 * An existing summary is not deleted on failure. This is not transactional.
	 * Constructor/start failures before handoff and arbitrary worker/protocol failures
	 * retain the same limits as makeChroms; an unstarted writer is never poisoned.
	 *
	 * @param genome Genome build number
	 * @param chroms Number of chromosomes to process
	 * @param name Genome name for metadata
	 * @param source Source file path for metadata
	 * @param unload Whether to unload chromosomes after processing
	 * @param scafNamePrefix Prefix flag recorded in summary metadata; does not rename scaffolds
	 * @throws RuntimeException If the info writer reports failure, before writing a new summary
	 * @deprecated Generate metadata during conversion rather than rereading serialized arrays
	 */
	@Deprecated
	public static void writeInfo(int genome, int chroms, String name, String source, boolean unload, boolean scafNamePrefix){
		Data.GENOME_BUILD=genome;
		Data.chromosomePlusMatrix=new ChromosomeArray[chroms+1];

		String outRoot=Data.ROOT_GENOME+genome+"/";
		TextStreamWriter info=null;
		boolean started=false;
		Throwable failure=null;
		long bases=0;
		long definedBases=0;
		long contigSum=0;

		//DNA-054: clean up a successfully started writer after producer failure.
		//TODO: Constructor/start failures before handoff may leave an inaccessible
		//or unstarted backend; the writer API has no close-without-start operation.
		try{
			info=new TextStreamWriter(outRoot+"info.txt", true, false, false);
			info.start();
			started=true;
			info.print("#Chromosome sizes\n");
			try{
				info.print("#Generated on\t"+new Date()+"\n");
			}catch(Exception e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			info.print("#Version\t"+VERSION+"\n");
			info.print("#chrom\tscaffolds\tcontigs\tlength\tdefined\tundefined\tstartPad\tstopPad\n");

			for(int chrom=1; chrom<=chroms; chrom++){
				ChromosomeArray ca=Data.getChromosome(chrom);
				int[] v=countInfo(ca);
				info.print(v[0]+"\t"+v[1]+"\t"+v[2]+"\t"+v[3]+"\t"+v[4]+"\t"+v[5]+"\t"+v[6]+"\t"+v[7]+"\n");

				bases+=v[3];
				definedBases+=v[4];
				contigSum+=v[2];
				if(unload){Data.unload(chrom, false);}
			}
		}catch(RuntimeException | Error problem){
			failure=problem;
		}finally{
			boolean poisoned=false;
			try{
				if(started){info.poison(); poisoned=true;}
			}catch(RuntimeException | Error problem){failure=retainFailure(failure, problem);}
			//Do not join when the termination request failed; arbitrary protocol/worker
			//failure recovery is outside this bounded producer-cleanup guarantee.
			try{
				if(poisoned){
					info.waitForFinish();
					if(info.errorState){
						failure=retainFailure(failure, new RuntimeException("Failed writing reference info: "+outRoot+"info.txt"));
					}
				}
			}catch(RuntimeException | Error problem){failure=retainFailure(failure, problem);}
		}
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		if(failure instanceof Error){throw (Error)failure;}

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
	}

	/**
	 * Core conversion method that processes FASTA file into ChromosomeArrays.
	 * Creates output directories, processes scaffolds into chromosomes with padding,
	 * generates metadata files, and optionally writes chromosome objects to disk.
	 * NODISK controls metadata destinations and directory cleanup, not writeChroms.
	 * With overwrite enabled, removes existing chromosome and selected index files
	 * before processing. This instance closes its constructed input and joins metadata
	 * writers whose start handshake returned, on normal and exceptional producer exits.
	 * The original producer exception is preserved; distinct cleanup failures are
	 * suppressed onto it. Retained errors reject before publishing a new summary.
	 * A throwing constructor/start may leave an inaccessible or unstarted backend;
	 * that startup limitation and arbitrary worker failures are not repaired here.
	 * Background chromosome writes retain ReadWrite's shared completion accounting
	 * when WAIT_FOR_WRITING is enabled, including after producer failure. This counter
	 * is not a failure report; IO042 and oversize KillSwitch behavior are unchanged.
	 *
	 * @param fname Input FASTA file path
	 * @param outRoot Output directory path ending in a separator
	 * @param genomeName Genome name for metadata
	 * @param genScaffolds Whether to generate scaffold information
	 * @param writeChroms Whether to write ChromosomeArray objects to disk
	 * @param r Mutable destination list, required when RETAIN=true
	 * @param scafNamePrefix Flag recorded in summary metadata; input header strings are unchanged
	 * @return Number of chromosomes created
	 * @throws RuntimeException If a retained input or metadata output error is detected
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

		ByteFile1 tf=null;
		int chrom=1;
		TextStreamWriter infoWriter=null, scafWriter=null;
		boolean infoStarted=false, scafStarted=false;
		boolean errorState=false;
		Throwable failure=null;
		ArrayList<String> infolist=null, scaflist=null;

		//DNA-024: only fully returned resources can be cleaned up through these APIs.
		//TODO: Constructor/start failures before handoff may leave an open backend;
		//TextStreamWriter has no public close-without-start operation. Never poison NEW.
		try{
			tf=new ByteFile1(fname, false);
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
				infoStarted=true;
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
					scafStarted=true;
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

		}catch(RuntimeException | Error problem){
			failure=problem;
		}finally{
			try{
				if(tf!=null){errorState|=tf.close();}
			}catch(RuntimeException | Error problem){failure=retainFailure(failure, problem);}
			boolean infoPoisoned=false, scafPoisoned=false;
			try{
				if(infoStarted){infoWriter.poison(); infoPoisoned=true;}
			}catch(RuntimeException | Error problem){failure=retainFailure(failure, problem);}
			try{
				if(scafStarted){scafWriter.poison(); scafPoisoned=true;}
			}catch(RuntimeException | Error problem){failure=retainFailure(failure, problem);}
			//Do not wait forever for a writer whose termination request itself failed.
			//Ordinary conversion uses no ordered jobs; arbitrary worker/protocol failures
			//remain outside this producer-cleanup repair.
			try{
				if(infoPoisoned){infoWriter.waitForFinish(); errorState|=infoWriter.errorState;}
			}catch(RuntimeException | Error problem){failure=retainFailure(failure, problem);}
			try{
				if(scafPoisoned){scafWriter.waitForFinish(); errorState|=scafWriter.errorState;}
			}catch(RuntimeException | Error problem){failure=retainFailure(failure, problem);}
			try{
				if(WAIT_FOR_WRITING && ReadWrite.countActiveThreads()>0){
					System.err.println("Waiting for writing to finish.");
					ReadWrite.waitForWritingToFinish();
					System.err.println("Finished.");
				}
			}catch(RuntimeException | Error problem){failure=retainFailure(failure, problem);}
		}
		if(errorState){
			failure=retainFailure(failure, new RuntimeException("Failed reading reference or writing reference metadata: "+fname+
					"; incomplete reference must not be used."));
		}
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		if(failure instanceof Error){throw (Error)failure;}

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

	/** Keeps the first producer/cleanup failure and attaches distinct later failures.
	 * @param first Existing unchecked failure, or null
	 * @param next Non-null unchecked failure from a subsequent cleanup operation
	 * @return The first failure if present, otherwise next */
	private static Throwable retainFailure(Throwable first, Throwable next){
		if(first==null){return next;}
		if(first!=next){first.addSuppressed(next);}
		return first;
	}

	/**
	 * Creates the next ChromosomeArray by reading scaffolds and applying padding.
	 * Merges multiple scaffolds into single chromosome if MERGE_SCAFFOLDS=true.
	 * Applies START_PADDING, MID_PADDING between scaffolds, and END_PADDING.
	 * MAX_LENGTH controls packing once a scaffold is present; a single larger scaffold
	 * may occupy its own chromosome up to MAX_SINGLE_SCAFFOLD. Both new and carried
	 * scaffolds are filtered by MIN_SCAFFOLD. A pending scaffold and its header remain
	 * on this instance when a packing boundary defers them to the next call.
	 * Counts/metadata use zero-based starts and base lengths; array maxIndex is inclusive.
	 *
	 * @param tf ByteFile1 for reading FASTA input
	 * @param chrom Chromosome number being created
	 * @param infoWriter TextStreamWriter for info.txt output
	 * @param scafWriter TextStreamWriter for scaffolds.txt.gz output
	 * @param infolist ArrayList for info data if NODISK=true
	 * @param scaflist Optional list to receive scaffold metadata in memory mode
	 * @return ChromosomeArray with scaffolds and padding, null if no more data
	 */
	private ChromosomeArray makeNextChrom(ByteFile1 tf, int chrom, TextStreamWriter infoWriter, TextStreamWriter scafWriter, ArrayList<String> infolist, ArrayList<String> scaflist){
		ChromosomeArray ca=new ChromosomeArray(chrom, (byte)Shared.PLUS, 0, 120000+START_PADDING);
		ca.maxIndex=-1;
		for(int i=0; i<START_PADDING; i++){ca.set(i, 'N');}

		if(verbose){System.err.println("chrom="+chrom+", lastHeader="+lastHeader+", nextHeader="+nextHeader);}

		int scaffolds=0;
		if(currentScaffold!=null && currentScaffold.length()>0){
			assert(lastHeader!=null);
			if(currentScaffold.length()>=MIN_SCAFFOLD){ //#002: honor minscaf for carried-over scaffolds too (was unconditional)
				assert(currentScaffold.length()+END_PADDING+ca.maxIndex<MAX_SINGLE_SCAFFOLD) : currentScaffold.length()+", "+END_PADDING+", "+ca.maxIndex+", "+MAX_SINGLE_SCAFFOLD;

//				System.err.println("A: Writing a scaffold because currentScaffold = "+currentScaffold);
				scaffoldSum++;
				if(scafWriter!=null){scafWriter.print(chrom+"\t"+scaffoldSum+"\t"+(ca.maxIndex+1)+"\t"+currentScaffold.length()+"\t"+lastHeader+"\n");}
				if(scaflist!=null && lastHeader!=null){
					scaflist.add(chrom+"\t"+scaffoldSum+"\t"+(ca.maxIndex+1)+"\t"+currentScaffold.length()+"\t"+lastHeader);
					if(verbose){System.err.println("A: Added to scaflist: "+scaflist.get(scaflist.size()-1));}
				}
				ca.set(ca.maxIndex+1, currentScaffold);
				scaffolds++;
			}

			currentScaffold.setLength(0);
			lastHeader=nextHeader;
		}

//		if()

		while((currentScaffold=nextScaffold(currentScaffold, tf))!=null){
			//#001: a single scaffold too large to fit even an empty chromosome — fail loudly instead of
			//silently building 0 chroms (first scaffold) or an over-MAX_LENGTH array (carried scaffold).
			//This existing KillSwitch policy predates producer finally cleanup: it prevented
			//non-daemon metadata writers from hanging after oversize rejection. Preserve the
			//whole-VM halt policy here; DNA-024 does not change oversize-input handling.
			if((long)currentScaffold.length()+START_PADDING+END_PADDING>MAX_SINGLE_SCAFFOLD){
				KillSwitch.kill("Scaffold '"+lastHeader+"' length "+currentScaffold.length()+
					" exceeds the maximum indexable single-chromosome length ~"+(MAX_SINGLE_SCAFFOLD-START_PADDING-END_PADDING)+
					" (2^31 int-array limit); cannot index this reference. Split the scaffold or use a smaller-genome workflow.");
			}
			//#003 (Eru 2026-06-25): the "scaffolds>0 &&" guard lets a single scaffold larger than MAX_LENGTH (but
			//<= MAX_SINGLE_SCAFFOLD) occupy its OWN chromosome instead of being broken-then-rejected.  No-op for any
			//scaffold that fits MAX_LENGTH (those never satisfy this on an empty chrom anyway), so normal multi-scaffold
			//packing -- and alignment output for ordinary genomes -- is unchanged.
			if(scaffolds>0 && currentScaffold.length()+MID_PADDING+END_PADDING+ca.maxIndex>MAX_LENGTH){break;}
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
			while(terminalN<=END_PADDING && ca.maxIndex<MAX_SINGLE_SCAFFOLD-1){
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

	/**
	 * Reads the next scaffold sequence from FASTA file.
	 * Accumulates sequence lines until next header or end of file.
	 * Updates nextHeader with the following header, excluding its leading greater-than
	 * character. Sequence lines are concatenated as bytes; normalization happens when
	 * the sequence is copied into a ChromosomeArray. This helper does not own input closure.
	 *
	 * @param sb ByteBuilder to reuse for sequence accumulation
	 * @param tf ByteFile1 for reading FASTA input
	 * @return Reused builder, possibly empty at a header boundary; null after exhaustion
	 */
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

	/** Header belonging to currentScaffold, without the leading greater-than character. */
	private String lastHeader;
	/** Header line for the next scaffold to be processed */
	private String nextHeader;
	/** Current sequence buffer, retained across calls when a scaffold is deferred by packing. */
	private ByteBuilder currentScaffold;
	/** Running count of total scaffolds processed across all chromosomes */
	private long scaffoldSum=0;
	/** Running sum of total bases across all chromosomes */
	private long lengthSum=0;
	/** Running sum of defined (non-N) bases across all chromosomes */
	private long definedSum=0;
	/** Running sum of undefined (N) bases across all chromosomes */
	private long undefinedSum=0;
	/** Running sum of total contigs across all chromosomes */
	private long contigSum=0;

	/** Returns the current version number for format compatibility */
	public static final int currentVersion(){return VERSION;}

	/** Whether to merge multiple scaffolds into single chromosomes */
	public static boolean MERGE_SCAFFOLDS=true;
	/** Retained legacy parsed state, restored when main2 returns or throws after preprocessing.
	 * Chromosome output calls the background writer regardless of this value. */
	public static boolean WRITE_IN_THREAD=false;
	/** Whether to overwrite existing output files */
	public static boolean overwrite=true;
	/** Legacy parsed option; not consulted by this implementation's output calls. */
	public static boolean append=false;
	/** Number of N bases to add at the beginning of each chromosome */
	public static int START_PADDING=8000; //Always applied
	/** Padding added before subsequent input scaffolds when a chromosome already contains a scaffold. */
	public static int MID_PADDING=300; //Applied when merging scaffolds
	/** Legacy terminal-padding threshold; extension uses an inclusive comparison and a length cap. */
	public static int END_PADDING=8000;
	/** Minimum scaffold length in bases to include in output */
	public static int MIN_SCAFFOLD=1;
	/** Minimum gap length in bases to trigger a new contig boundary */
	public static int contigTrigger=10;
	/** Format version number for compatibility tracking */
	public static int VERSION=5;
	/** Multi-scaffold packing threshold; a permitted single scaffold may exceed this value. */
	public static int MAX_LENGTH=(1<<29)-200000;
	/** Maximum allowed length for a single (unmergeable) scaffold, in bases.  Defaults to MAX_LENGTH so behavior
	 * is byte-for-byte unchanged for every existing caller.  BBMap5 raises it (to ~MAX_ARRAY_LEN/2) so a
	 * wheat-scale chromosome becomes its own chrom rather than being rejected, while MAX_LENGTH still governs
	 * multi-scaffold packing -- so packed-chrom layout, and thus alignment output, stays identical for normal
	 * genomes.  Only takes effect when set larger than MAX_LENGTH. */
	public static int MAX_SINGLE_SCAFFOLD=MAX_LENGTH;
//	public static boolean TRANSLATE_U_TO_T;

	/** Whether to print verbose debugging information */
	public static boolean verbose=false;
	/** Whether to retain ChromosomeArray objects in memory for return */
	public static boolean RETAIN=false;
	/** Waits for the shared background object-writer counter; main2 also enables it for gzip/pigz. */
	public static boolean WAIT_FOR_WRITING=true;
	/** Redirects metadata to memory and skips directory cleanup; writeChroms remains independent. */
	public static boolean NODISK=false;
	/** Copied from Data.GENOME_BUILD when publishing lists; may differ from main2's build argument. */
	public static int LISTBUILD=-1;
	/** Shared info, scaffold and summary lines for memory-mode conversion, updated on normal completion. */
	public static ArrayList<String> INFO_LIST, SCAF_LIST, SUMMARY_LIST;

//	public static boolean GENERATE_SCAFFOLD_INFO=true;

}
