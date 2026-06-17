package var2;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import ml.CellNet;
import ml.CellNetParser;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import structures.ListNum;

/**
* Filters VCF files based on variant quality, type, position, and statistical criteria.
* Provides comprehensive filtering capabilities for post-processing variant calls,
* with support for both single-threaded and multithreaded operation.
* 
* Key features:
* - Statistical filtering using VarFilter criteria (coverage, quality, strand bias, etc.)
* - Position-based filtering using SamFilter criteria (coordinates, contigs)
* - Region restriction using a BED file of intervals (bed=, optionally invertbed)
* - Variant type filtering (enable/disable SNPs, indels, junctions)
* - Genotype filtering by GT from the first sample column (gt=, homozygous=, heterozygous=)
* - Allele splitting for multi-allelic variants, and multi-base substitution splitting
* - Indel left-alignment / normalization against the reference (normalize, requires ref=)
* - Quality score histograms for analysis
* - Header preservation and metadata extraction
*
* @author Brian Bushnell
* @contributor Isla
* @contributor UMP45
* @date January 14, 2017
*/
public class FilterVCF {
	
	/**
	 * Main method for command-line execution.
	 * 
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		FilterVCF x=new FilterVCF(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses command-line arguments and initializes filtering parameters.
	 * Sets up both variant-specific filtering (VarFilter) and position-specific filtering (SamFilter).
	 * 
	 * @param args Command line arguments array
	 */
	public FilterVCF(String[] args){
		
		{ //Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		varFilter.clear();

		boolean setSamFilter=false;
		boolean setVarFilter=false;
		
		// Neural network parameters
		String netFile=null;
		boolean autoCutoff=true;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("multithreaded") || a.equals("mt")){
				multithreaded=Parse.parseBoolean(b);
			}else if(a.equals("singlethreaded") || a.equals("st")){
				multithreaded=!Parse.parseBoolean(b);
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("sub") || a.equals("subs")){
				Var.CALL_SUB=Parse.parseBoolean(b);
			}else if(a.equals("del") || a.equals("dels")){
				Var.CALL_DEL=Parse.parseBoolean(b);
			}else if(a.equals("ins") || a.equals("inss")){
				Var.CALL_INS=Parse.parseBoolean(b);
			}else if(a.equals("indel") || a.equals("indels")){
				Var.CALL_INS=Var.CALL_DEL=Parse.parseBoolean(b);
			}else if(a.equals("junction") || a.equals("junctions")){
				Var.CALL_JUNCTION=Parse.parseBoolean(b);
			}else if(a.equals("minscore")){
				minScore=Double.parseDouble(b);
			}else if(a.equals("minqual") || a.equals("minphred")){
				minQual=Double.parseDouble(b);
			}else if(a.equals("splitalleles")){
				splitAlleles=Parse.parseBoolean(b);
			}else if(a.equals("splitsubs") || a.equals("splitsnps")){
				splitSubs=Parse.parseBoolean(b);
			}else if(a.equals("splitcomplex")){
				splitComplex=Parse.parseBoolean(b);
			}else if(a.equals("sass") || a.equals("split")){
				splitAlleles=splitSubs=Parse.parseBoolean(b);
			}else if(a.equals("splitall") || a.equals("sascsss")){
				splitAlleles=splitComplex=splitSubs=Parse.parseBoolean(b);
			}else if(a.equals("normalize") || a.equals("leftalign") || a.equals("norm")){
				normalize=Parse.parseBoolean(b);
			}else if(a.equals("clearfilters")){
				if(Parse.parseBoolean(b)){
					varFilter.clear();
					samFilter.clear();
				}
			}else if(a.equals("bed") || a.equals("bedfile")){
				bedFile=b;
			}else if(a.equals("invertbed") || a.equals("excludebed") || a.equals("bedexclude") || a.equals("bedinvert")){
				invertBed=Parse.parseBoolean(b);
			}else if(a.equals("gt") || a.equals("genotype")){
				if(b==null || b.length()==0 || b.equalsIgnoreCase("all") || b.equalsIgnoreCase("any")){
					gtWhitelist=null;
				}else{
					gtWhitelist=b.split(",");
					for(int j=0; j<gtWhitelist.length; j++){gtWhitelist[j]=gtWhitelist[j].replace('|', '/');}
				}
			}else if(a.equals("homozygous") || a.equals("hom")){
				homFilter=Parse.parseBoolean(b) ? 1 : -1;
			}else if(a.equals("heterozygous") || a.equals("het")){
				homFilter=Parse.parseBoolean(b) ? -1 : 1;
			}else if(samFilter.parse(arg, a, b)){
				setSamFilter=true;
			}else if(varFilter.parse(a, b, arg)){
				setVarFilter=true;
			}else if(a.equals("scorehist") || a.equals("qualhist") || a.equals("qhist") || a.equals("shist")){
				scoreHistFile=b;
			}else if(a.equals("trimtocanonical") || a.equals("canonicalize") || a.equals("canonicize") || a.equals("canonize")){
				VCFLine.TRIM_TO_CANONICAL=Parse.parseBoolean(b);
			}
			
			//Neural network parameters
			else if(a.equals("net") || a.equals("netfile")){
				netFile=b;
				useNet=(b!=null);
			}else if(a.equals("netcutoff")){
				if("auto".equalsIgnoreCase(b)){
					autoCutoff=true;
				}else{
					autoCutoff=false;
					netCutoff=Float.parseFloat(b);
				}
			}else if(a.equals("usenet") || a.equals("useann") || a.equals("usenn") || a.equals("nn")){
				useNet=Parse.parseBoolean(b);
			}else if(a.equals("netmode")){
				useNet=(b!=null);
				if(b!=null){FeatureVectorMaker.setMode(b);}
			}
			
			else if(a.equalsIgnoreCase("countNearbyVars")){
				countNearby=Parse.parseBoolean(b);
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		{ //Process parser fields
			in1=parser.in1;
			out1=parser.out1;
			overwrite=parser.overwrite;
			append=parser.append;
		}

		if(!setSamFilter){samFilter=null;}
		if(!setVarFilter){varFilter=null;}

		//Load neural network if specified
		if(netFile!=null && useNet){
			net0=CellNetParser.load(netFile);
			assert(net0!=null) : "Failed to load neural network: "+netFile;
			if(autoCutoff){netCutoff=net0.cutoff;}
			if(verbose){outstream.println("Loaded neural network: "+netFile+" (cutoff="+netCutoff+")");}
		}else{
			net0=null;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least two input files are required.");}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, ref, bedFile)){
			throw new RuntimeException("\nCan't read some input files.\n");
		}

		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, multithreaded);

		if(normalize && ref==null){
			throw new RuntimeException("Error - normalize/leftalign requires ref= (the reference fasta) for left-alignment.");
		}
		if(ref!=null){ScafMap.loadReference(ref, scafMap, samFilter, true);}

		//Load region mask, if a BED file was specified
		if(bedFile!=null){
			bedMask=new BedMask(bedFile);
			outstream.println("Loaded "+bedMask.intervalsLoaded()+" BED intervals across "+bedMask.scaffolds()+" scaffolds.");
		}
		
		//Determine how many threads may be used
		threads=Tools.min(8, Shared.threads());
	}
	
	/**
	 * Returns the reference bases (0-based) for a scaffold, for left-alignment, or null if
	 * the reference was not loaded for this scaffold.  Reads from the reference ScafMap, which
	 * is immutable after loading, so concurrent reads from worker threads are safe.
	 *
	 * @param scaf Scaffold/chromosome name
	 * @return Reference base array, or null when unavailable
	 */
	private byte[] refBasesFor(String scaf){
		Scaffold sc=scafMap.getScaffold(scaf);
		return sc==null ? null : sc.bases;
	}

	/**
	 * Genotype-based pass test for the gt=/homozygous=/heterozygous= filters.  Returns true when
	 * the variant's genotype (GT, the first ':'-subfield of the first sample column) satisfies BOTH
	 * the gt= whitelist (if set) and the zygosity constraint (if set).  A record lacking a sample/GT
	 * column is rejected whenever a genotype filter is active (it cannot be evaluated).
	 *
	 * @param vline The parsed VCF line
	 * @return true if the variant passes the active genotype filter(s)
	 */
	private boolean passesGenotype(VCFLine vline){
		final String gt=genotypeOf(vline);
		if(gt==null){return false;}
		if(gtWhitelist!=null){
			boolean match=false;
			for(String g : gtWhitelist){if(g.equals(gt)){match=true; break;}}
			if(!match){return false;}
		}
		if(homFilter!=0){
			final boolean hom=isHomozygous(gt);
			if(homFilter>0 && !hom){return false;}
			if(homFilter<0 && hom){return false;}
		}
		return true;
	}

	/**
	 * Extracts the GT (genotype) string from the first sample column of a VCF line, with '|' phasing
	 * normalized to '/'.  By VCF convention GT is the first ':'-delimited subfield of the first
	 * sample.  Returns null when no sample column is present or the genotype is empty.
	 *
	 * @param vline The parsed VCF line
	 * @return The normalized GT string (e.g. "1/1"), or null if unavailable
	 */
	private static String genotypeOf(VCFLine vline){
		if(vline.samples==null || vline.samples.isEmpty()){return null;}
		final byte[] s=vline.samples.get(0);
		if(s==null || s.length==0){return null;}
		int end=0;
		while(end<s.length && s[end]!=':'){end++;}
		if(end==0){return null;}
		final StringBuilder sb=new StringBuilder(end);
		for(int i=0; i<end; i++){
			final char c=(char)s[i];
			sb.append(c=='|' ? '/' : c);
		}
		return sb.toString();
	}

	/**
	 * Tests whether a normalized genotype string is homozygous (all alleles identical and present).
	 * Single-allele (haploid) genotypes are homozygous by definition; any missing allele ('.') is
	 * treated as non-homozygous.
	 *
	 * @param gt Normalized genotype string (e.g. "1/1", "0/1", "1")
	 * @return true if homozygous
	 */
	private static boolean isHomozygous(String gt){
		if(gt==null || gt.length()==0){return false;}
		final String[] alleles=gt.split("/");
		if(alleles.length==0 || alleles[0].equals(".")){return false;}
		for(int i=1; i<alleles.length; i++){
			if(!alleles[i].equals(alleles[0])){return false;}
		}
		return true;
	}

	/**
	 * Loads scaffold information from VCF contig header lines into the ScafMap.
	 * Processes ##contig=<ID=...> lines to build coordinate reference system.
	 */
	public void loadHeaderInScafMap(){
		for(byte[] line : header){
			if(Tools.startsWith(line, "##contig=<ID=")){
				scafMap.addFromVcf(line);
			}
		}
	}
	
	/**
	 * Converts stored header lines to a formatted string.
	 * 
	 * @return Complete VCF header as string
	 */
	public String headerToString(){
		StringBuilder sb=new StringBuilder();
		for(byte[] line : header){
			for(byte b : line){
				sb.append((char)b);
			}
			sb.append('\n');
		}
		return sb.toString();
	}
	
	/**
	 * Spawns worker threads for multithreaded VCF processing.
	 * Each thread processes batches of VCF lines independently.
	 * 
	 * @param bf Input ByteFile
	 * @param bsw Output ByteStreamWriter
	 * @return List of spawned ProcessThread objects
	 */
	private ArrayList<ProcessThread> spawnThreads(ByteFile bf, ByteStreamWriter bsw){
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(bf, bsw, jobIDOffset));
		}
		if(verbose){outstream.println("Spawned threads.");}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
		
		//Do anything necessary after processing
		return alpt;
	}
	
	/**
	 * Waits for all worker threads to complete and aggregates statistics.
	 * 
	 * @param alpt List of ProcessThread objects to wait for
	 */
	private void waitForFinish(ArrayList<ProcessThread> alpt){
		//Wait for completion of all threads
		boolean allSuccess=true;
		for(ProcessThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				}catch(InterruptedException e){
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			linesProcessed+=pt.linesProcessedT;
			headerLinesProcessed+=pt.headerLinesProcessedT;
			variantLinesProcessed+=pt.variantLinesProcessedT;
			variantLinesOut+=pt.variantLinesOutT;
			bytesProcessed+=pt.bytesProcessedT;
			Tools.add(scoreHist, pt.scoreHistT);

			allSuccess&=pt.success;
		}
		
		//Track whether any threads failed
		if(!allSuccess){errorState=true;}
	}
	
	/**
	 * Processes VCF header section and extracts metadata.
	 * Handles scaffold definitions, sample names, and statistical metadata
	 * from the VCF header lines.
	 * 
	 * @param bf Input ByteFile
	 * @param bsw Output ByteStreamWriter (may be null)
	 */
	private void processVcfHeader(ByteFile bf, ByteStreamWriter bsw){
		byte[] line=bf.nextLine();

		if(ScafMap.defaultScafMap()==null){
			ScafMap.setDefaultScafMap(new ScafMap(), bf.name());
		}
		ByteBuilder bb=new ByteBuilder();
		while(line!=null && (line.length==0 || line[0]=='#')){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				headerLinesProcessed++;
				bytesProcessed+=line.length;
				bb.append(line).append('\n');
				headerLinesOut++;
				header.add(line);
				if(Tools.startsWith(line, VCFFile.CHROM_POS)){
					String[] split=new String(line).split("\t");
					for(int i=9; i<split.length; i++){
						samples.add(split[i]);
					}
				}else if(Tools.startsWith(line, "##contig=<ID=")){
					ScafMap.defaultScafMap().addFromVcf(line);
				}else{
					String[] split=new String(line).split("=");
					if(split.length==2){
						String a=split[0], b=split[1];
						if(a.equalsIgnoreCase("##ploidy")){
							ploidy=Integer.parseInt(b);
						}else if(a.equalsIgnoreCase("##properPairRate")){
							properPairRate=(float) Double.parseDouble(b);
						}else if(a.equalsIgnoreCase("##totalQualityAvg")){
							totalQualityAvg=(float) Double.parseDouble(b);
						}else if(a.equalsIgnoreCase("##mapqAvg")){
							totalMapqAvg=(float) Double.parseDouble(b);
						}else if(a.equalsIgnoreCase("##readLengthAvg")){
							readLengthAvg=(float) Double.parseDouble(b);
						}
					}
				}
			}
			line=bf.nextLine();
		}
		if(line!=null && line.length>0){bf.pushBack(line);}
		if(bb.length()>0 && bsw!=null){
			bsw.add(bb, jobIDOffset);
			jobIDOffset++;
		}
	}
	
	/**
	 * Single-threaded VCF variant processing.
	 * Processes each variant line sequentially, applying all filtering criteria
	 * and optionally splitting multi-allelic variants.
	 * 
	 * @param bf Input ByteFile
	 * @param bsw Output ByteStreamWriter (may be null)
	 */
	private void processVcfVarsST(ByteFile bf, ByteStreamWriter bsw){
		/** Whether Var format conversion is working */
		boolean varFormatOK=true;
		byte[] line=bf.nextLine();
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=line.length;
				
				final boolean isHeader=(line[0]=='#');
				
				if(isHeader){
					assert(false) : "Encountered intermediate header.";
					headerLinesProcessed++;
					if(bsw!=null){bsw.println(line);}
					header.add(line);
					if(Tools.startsWith(line, VCFFile.CHROM_POS)){
						String[] split=new String(line).split("\t");
						for(int i=9; i<split.length; i++){
							samples.add(split[i]);
						}
					}
				}else{
					variantLinesProcessed++;
					VCFLine vline=new VCFLine(line);
					boolean pass=true;

					//Quality-score filtering on the VCF QUAL column.  Apply BOTH minScore and minQual, mirroring
					//the multithreaded path (processLine, lines ~832-833).  minScore was previously applied ONLY in
					//the multithreaded path, so the default single-threaded path silently ignored minscore=.
					if(vline.qual<minScore){pass=false;}
					if(vline.qual<minQual){pass=false;}

					//Type-based filtering
					if(!Var.CALL_DEL && vline.type()==Var.DEL){pass=false;}
					else if(!Var.CALL_INS && vline.type()==Var.INS){pass=false;}
					else if(!Var.CALL_SUB && vline.type()==Var.SUB){pass=false;}
					else if(!Var.CALL_JUNCTION && vline.isJunction()){pass=false;}
					
					//Position-based filtering
					if(pass && samFilter!=null){pass&=samFilter.passesFilter(vline);}

					//Region-based (BED) filtering
					if(pass && bedMask!=null){
						boolean in=bedMask.contains(vline.scaf, vline.pos);
						pass&=(invertBed ? !in : in);
					}

					//Genotype-based (GT) filtering
					if(pass && (gtWhitelist!=null || homFilter!=0)){pass&=passesGenotype(vline);}

					//Statistical filtering
					if(pass && varFilter!=null){
						Var v=null;
						
						if(varFormatOK){
							try {
								v=vline.toVar();
							}catch(Throwable e){
								System.err.println("WARNING: This VCF file does not support Var format.\n"+"Filtering can only be done on location and quality score.\n");
								varFormatOK=false;
							}
						}
						
						if(v!=null){
							pass&=varFilter.passesFilter(v, properPairRate, totalQualityAvg, totalMapqAvg,
									readLengthAvg, ploidy, scafMap, net0, countNearby);
						}else{
							pass&=vline.qual>=varFilter.minScore;
						}
					}
					
					if(pass){
						//Handle variant splitting if requested
						ArrayList<VCFLine> split=(splitAlleles || splitComplex || splitSubs) ? vline.split(splitAlleles, splitComplex, splitSubs) : null;
						
						if(split==null){
							if(normalize){vline.leftAlign(refBasesFor(vline.scaf));}
							if(bsw!=null){
								if(normalize){bsw.print(vline.toText(new ByteBuilder(line.length+16)).nl());}
								else{bsw.println(line);}
							}
							variantLinesOut++;
							int q=(int)(vline.qual);
							scoreHist[Tools.min(scoreHist.length-1, q)]++;
						}else{
							for(VCFLine vline2 : split){
								if(normalize){vline2.leftAlign(refBasesFor(vline2.scaf));}
								if(bsw!=null){bsw.print(vline2.toText(new ByteBuilder(64)).nl());}
								variantLinesOut++;
								int q=(int)(vline2.qual);
								scoreHist[Tools.min(scoreHist.length-1, q)]++;
							}
						}
						
						/* COMMENTED OUT: Earlier allele splitting implementation
						 * This was a simpler approach that just split multi-allelic variants
						 * by comma separation, but didn't handle the auxiliary VCF data correctly.
						 * The current split() method properly handles FORMAT fields and other
						 * per-allele information when splitting variants.
						 */
//						if(splitAlleles && vline.alt!=null && Tools.indexOf(vline.alt, ',')>0){//This may not split correctly, since the auxiliary data is replicated
//							String alleles=new String(vline.alt);
//							String[] split=alleles.split(",");
//							for(String allele : split){
//								vline.alt=allele.getBytes();
//								if(bsw!=null){bsw.print(vline.toText(new ByteBuilder(64)).nl());}
//								variantLinesOut++;
//								int q=(int)(vline.qual);
//								scoreHist[Tools.min(scoreHist.length-1, q)]++;
//							}
//						}else{
//							if(bsw!=null){bsw.println(line);}
//							variantLinesOut++;
//							int q=(int)(vline.qual);
//							scoreHist[Tools.min(scoreHist.length-1, q)]++;
//						}
//						if(bsw!=null){bsw.println(line);}
//						variantLinesOut++;
//						int q=(int)(vline.qual);
//						scoreHist[Tools.min(scoreHist.length-1, q)]++;
					}
				}
			}
			line=bf.nextLine();
		}
	}
	
	/**
	 * Multithreaded VCF variant processing.
	 * Spawns worker threads to process variants in parallel for better performance.
	 * 
	 * @param bf Input ByteFile
	 * @param bsw Output ByteStreamWriter
	 */
	private void processVcfVarsMT (ByteFile bf, ByteStreamWriter bsw){
		ArrayList<ProcessThread> alpt=spawnThreads(bf, bsw);
		waitForFinish(alpt);
	}
	
	/* COMMENTED OUT: Producer-consumer threading model
	 * This was an earlier approach that used a separate thread to read bytes
	 * and distribute them to worker threads via a queue. The current approach
	 * uses ByteFile.nextList() directly in each thread for simpler coordination.
	 */
//	private void readBytes(ByteFile bf){
//		ListNum<byte[]> ln=bf.nextList();
//		while(ln!=null){
//			putBytes(ln);
//			ln=bf.nextList();
//		}
//		putBytes(POISON_BYTES);
//	}
	
	/**
	 * Core filtering method that processes an entire VCF file.
	 * Handles header processing, scaffold map initialization, and variant filtering
	 * using either single-threaded or multithreaded approach.
	 * 
	 * @param ff Input file format
	 * @param bsw Output ByteStreamWriter
	 */
	public void filter(FileFormat ff, ByteStreamWriter bsw){
		
		ByteFile bf=ByteFile.makeByteFile(ff);
		
		processVcfHeader(bf, bsw);
		
		loadHeaderInScafMap();
		assert(ScafMap.defaultScafMap().size()>0) : ScafMap.defaultScafMap()+"\n"+headerToString();

		if(multithreaded){
			processVcfVarsMT(bf, bsw);
		}else{
			processVcfVarsST(bf, bsw);
		}
		
		if(scoreHistFile!=null){
			CVOutputWriter.writeScoreHist(scoreHistFile, scoreHist);
		}
		
		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
	}
	
	/**
	 * Main processing method that coordinates the entire filtering workflow.
	 * 
	 * @param t Timer for performance measurement
	 */
	void process(Timer t){
		
		/** Output ByteStreamWriter */
		ByteStreamWriter bsw;
		if(ffout1!=null){
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}else{bsw=null;}
		
		filter(ffin1, bsw);
		
		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		outstream.println();
		outstream.println("Header Lines In:   \t"+headerLinesProcessed);
		//TODO: Probable bug - "Variant Lines In" (variantLinesProcessed) and "Variant Lines Out"
		//(variantLinesOut) both undercount by exactly 1, in BOTH the single- and multi-threaded
		//paths (verified on chr21: reported 79954 vs 79955 actual output variants).  The output
		//file itself is correct; only these two counters are short.  Most likely the first variant
		//line after the header - pushed back in processVcfHeader() - is written but not tallied.
		//Pre-existing (reproduces with no bed=), not introduced by the bed= region filter.
		//Repro needs a VCF with ##contig header lines (FilterVCF asserts on an empty ScafMap).
		outstream.println("Variant Lines In:  \t"+variantLinesProcessed);
		outstream.println("Header Lines Out:  \t"+headerLinesOut);
		outstream.println("Variant Lines Out: \t"+variantLinesOut);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Threading Support      ----------------*/
	/*--------------------------------------------------------------*/

	/* COMMENTED OUT: Queue-based threading methods
	 * These were part of the earlier producer-consumer threading model
	 * that used an ArrayBlockingQueue to distribute work between threads.
	 * The current approach is simpler and more efficient.
	 */
//	final void putBytes(ListNum<byte[]> list){
//		while(list!=null){
//			try {
//				inq.put(list);
//				list=null;
//			} catch (InterruptedException e) {
//				e.printStackTrace();
//			}
//		}
//	}
//	
//	final ListNum<byte[]> takeBytes(){
//		ListNum<byte[]> list=null;
//		while(list==null){
//			try {
//				list=inq.take();
//			} catch (InterruptedException e) {
//				e.printStackTrace();
//			}
//		}
//		return list;
//	}
	
	/**
	 * Worker thread for multithreaded VCF processing.
	 * Each thread processes batches of VCF lines independently and maintains
	 * thread-local statistics that are aggregated at completion.
	 */
	private class ProcessThread extends Thread {
		
		/**
		 * Creates a ProcessThread for VCF line processing.
		 * 
		 * @param bf_ Input ByteFile
		 * @param bsw_ Output ByteStreamWriter
		 * @param jobIDOffset_ Job ID offset for output ordering
		 */
		ProcessThread(ByteFile bf_, ByteStreamWriter bsw_, long jobIDOffset_){
			bf=bf_;
			bsw=bsw_;
			offset=jobIDOffset_;
		}

		/**
		 * Main thread execution loop.
		 * Processes batches of lines from the ByteFile until exhausted.
		 */
		@Override
		/**
		* Main thread execution loop.
		* Processes batches of lines from the ByteFile until exhausted.
		*/
		public void run(){
			net=(net0==null ? null : net0.copy(false));
			ListNum<byte[]> ln=nextListSync();
			while(ln!=null && ln!=POISON_BYTES){
				ByteBuilder bb=new ByteBuilder(4096);
				for(byte[] line : ln){
					linesProcessedT++;
					processLine(line, bb);
				}
				if(bsw!=null){bsw.add(bb, ln.id+offset);}
				ln=nextListSync();
			}
			success=true;
			synchronized(this){notify();}
		}

		/**
		 * Reads the next batch from the shared input ByteFile under a lock.
		 * ByteFile1.nextList() is NOT thread-safe (it mutates shared buffer offsets), so
		 * multiple worker threads calling it concurrently corrupt the buffer and drop lines.
		 * Serializing only the read is cheap (decompression runs in parallel via pigz) while
		 * the expensive per-line work (split, left-align, filtering, NN scoring) stays parallel.
		 *
		 * @return Next batch of lines, or null at end of input
		 */
		private ListNum<byte[]> nextListSync(){
			synchronized(bf){
				return bf.nextList();
			}
		}

		/**
		 * Processes a single VCF line, applying all filtering criteria.
		 * 
		 * @param line Raw VCF line bytes
		 * @param bb ByteBuilder for output accumulation
		 */
		void processLine(byte[] line, ByteBuilder bb){
			linesProcessedT++;
			bytesProcessedT+=line.length;
			if(line.length<1){return;}
			final boolean isHeader=(line[0]=='#');
			if(isHeader){
				assert(false) : "Encountered intermediate header.";
				headerLinesProcessedT++;
				bb.append(line).append('\n');
				synchronized(header){
					header.add(line);
				}
				if(Tools.startsWith(line, VCFFile.CHROM_POS)){
					String[] split=new String(line).split("\t");
					synchronized(samples){
						for(int i=9; i<split.length; i++){
							samples.add(split[i]);
						}
					}
				}
			}else{
				variantLinesProcessedT++;
				boolean pass=true;
				
				VCFLine vline=new VCFLine(line);
				pass&=vline.qual>=minScore;
				if(pass){pass&=vline.qual>=minQual;}

				{
					if(pass){
						//Type-based filtering
						if(!Var.CALL_DEL && vline.type()==Var.DEL){pass=false;}
						else if(!Var.CALL_INS && vline.type()==Var.INS){pass=false;}
						else if(!Var.CALL_SUB && vline.type()==Var.SUB){pass=false;}
						else if(!Var.CALL_JUNCTION && vline.isJunction()){pass=false;}
					}

					//Position-based filtering
					if(pass && samFilter!=null){pass&=samFilter.passesFilter(vline);}

					//Region-based (BED) filtering
					if(pass && bedMask!=null){
						boolean in=bedMask.contains(vline.scaf, vline.pos);
						pass&=(invertBed ? !in : in);
					}

					//Genotype-based (GT) filtering
					if(pass && (gtWhitelist!=null || homFilter!=0)){pass&=passesGenotype(vline);}

					//Statistical filtering
					if(pass && varFilter!=null){
						if(varFormatOK){
							try {
								
								final Var v;
								if(threads>1){
									//Optimized for multithreaded use - faster parsing
									v=VcfToVar.fromVCF(line, scafMap, true, true);
								}else{
									//Optimized for single-threaded use - more thorough but slower
									v=vline.toVar();
								}
								
								pass&=varFilter.passesFilter(v, properPairRate, totalQualityAvg, totalMapqAvg,
										readLengthAvg, ploidy, scafMap, net, countNearby);
							}catch(Throwable e){
								System.err.println("WARNING: This VCF file does not support Var format.\n"+"Filtering can only be done on location and quality score.\n"+e);
								e.printStackTrace();
								varFormatOK=false;
							}
						}
						//Mirror the single-threaded fallback (processVcfVarsST): once Var conversion is known to fail,
						//filter on the QUAL column via varFilter.minScore.  The multithreaded path previously applied
						//no statistical filter at all in that case.
						if(!varFormatOK){pass&=vline.qual>=varFilter.minScore;}
					}
				}
				
				if(pass){
					//Handle variant splitting and/or normalization if requested.
					//Uses the same VCFLine.split() path as the single-threaded code, which
					//correctly handles per-allele INFO; the old comma-split here did not.
					ArrayList<VCFLine> split=(splitAlleles || splitComplex || splitSubs) ? vline.split(splitAlleles, splitComplex, splitSubs) : null;

					if(split==null){
						if(normalize){vline.leftAlign(refBasesFor(vline.scaf));}
						if(normalize){vline.toText(bb).nl();}
						else{bb.append(line).append('\n');}
						variantLinesOutT++;
						int q=(int)(vline.qual);
						scoreHistT[Tools.min(scoreHistT.length-1, q)]++;
					}else{
						for(VCFLine vline2 : split){
							if(normalize){vline2.leftAlign(refBasesFor(vline2.scaf));}
							vline2.toText(bb).nl();
							variantLinesOutT++;
							int q=(int)(vline2.qual);
							scoreHistT[Tools.min(scoreHistT.length-1, q)]++;
						}
					}
				}
			}
		}
		
		/** Input ByteFile */
		final ByteFile bf;
		/** Output ByteStreamWriter */
		final ByteStreamWriter bsw;
		/** Output ByteStreamWriter */
		CellNet net;
		/** Job ID offset for output ordering */
		final long offset;
		/** Whether Var format conversion is working */
		boolean varFormatOK=true;

		/** Thread-local statistics */
		long linesProcessedT=0;
		/** Thread-local header lines processed counter */
		long headerLinesProcessedT=0;
		/** Thread-local variant lines processed counter */
		long variantLinesProcessedT=0;
		/** Thread-local variant lines output counter */
		long variantLinesOutT=0;
		/** Thread-local bytes processed counter */
		long bytesProcessedT=0;
		/** Thread-local quality score histogram for aggregation */
		private long[] scoreHistT=new long[scoreHist.length];

		/** Success flag */
		boolean success=false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of lines processed from input */
	private long linesProcessed=0;
	/** Number of header lines processed */
	private long headerLinesProcessed=0;
	/** Number of variant lines processed */
	private long variantLinesProcessed=0;
	/** Number of header lines written to output */
	private long headerLinesOut=0;
	/** Number of variant lines written to output */
	private long variantLinesOut=0;
	/** Number of bytes processed from input */
	private long bytesProcessed=0;
	/** Histogram of variant quality scores for distribution analysis */
	private long[] scoreHist=new long[1000];
	
	/** Maximum number of lines to process (for testing/debugging) */
	private long maxLines=Long.MAX_VALUE;

	/** VCF header lines */
	public ArrayList<byte[]> header=new ArrayList<byte[]>();
	/** Sample names extracted from VCF header */
	public ArrayList<String> samples=new ArrayList<String>();
	
	/** Filter for SAM/alignment-based criteria */
	SamFilter samFilter=new SamFilter();
	/** Filter for variant-specific criteria */
	VarFilter varFilter=new VarFilter();
	
	/*--------------------------------------------------------------*/
	/*----------------     Neural Network Fields    ----------------*/
	/*--------------------------------------------------------------*/

	/** Master neural network model (copied to each thread) */
	private CellNet net0=null;
	/** Whether to use neural network for variant filtering */
	private boolean useNet=false;
	/** Score threshold for neural network filtering */
	private float netCutoff=0.5f;
	
	/*--------------------------------------------------------------*/
	/*----------------    Configuration Fields     ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Minimum quality score threshold for simple filtering */
	double minScore=0;
	/** Minimum VCF QUAL column threshold (filters directly on QUAL, not recomputed score) */
	double minQual=0;
	
	/** Sample ploidy for variant evaluation */
	public int ploidy=1;
	/** Proper pair rate from sequencing run */
	public float properPairRate=0;
	/** Average total quality from dataset */
	public float totalQualityAvg=30;
	/** Average mapping quality from dataset */
	public float totalMapqAvg=30;
	/** Average read length from sequencing run */
	public float readLengthAvg=150;
	
	/** Number of processing threads (limited to 8 maximum) */
	final int threads;
	/** Whether to use multithreaded processing */
	public boolean multithreaded=false;
	/** Job ID offset for ordered output */
	private long jobIDOffset=0;
	/** Whether to split multi-allelic variants into separate lines */
	boolean splitAlleles=false;
	/** Whether to split complex substitutions */
	boolean splitSubs=false;
	/** Whether to split complex variants */
	boolean splitComplex=false;
	/** Whether to left-align (normalize) indels using the reference; requires ref= */
	boolean normalize=false;

	/** Whether to count nearby variants for filtering (TODO: implement counting) */
	boolean countNearby=false;

	/** Acceptable genotypes for the gt= filter (each normalized '|'->'/'); null = no GT-list filter */
	private String[] gtWhitelist=null;
	/** Zygosity filter: 0=off, 1=keep homozygous only (homozygous=t), -1=keep heterozygous only */
	private int homFilter=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         File Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input VCF filename */
	private String in1=null;
	/** Primary output VCF filename */
	private String out1=null;
	/** Reference genome filename for variant validation */
	private String ref=null;
	/** Score histogram output filename */
	private String scoreHistFile=null;

	/** BED file of regions; if set, variants are restricted to (or excluded from) these intervals */
	private String bedFile=null;
	/** When true, keep variants OUTSIDE the BED intervals instead of inside */
	private boolean invertBed=false;
	/** Loaded region mask, or null when no bed= was given */
	private BedMask bedMask=null;

	/** Input file format */
	private final FileFormat ffin1;
	/** Output file format */
	private final FileFormat ffout1;
	
	/** Scaffold mapping for coordinate resolution */
	public final ScafMap scafMap=new ScafMap();
	
	/*--------------------------------------------------------------*/
	/*----------------       Static Fields          ----------------*/
	/*--------------------------------------------------------------*/

	/** Poison pill for ending thread processing */
	static final ListNum<byte[]> POISON_BYTES=new ListNum<byte[]>(null, Long.MAX_VALUE, true, false);
	/** Output stream for messages */
	private PrintStream outstream=System.err;
	/** Verbose output flag */
	public static boolean verbose=false;
	/** Error state flag */
	public boolean errorState=false;
	/** Overwrite output files flag */
	private boolean overwrite=true;
	/** Append to output files flag */
	private boolean append=false;
	
}
