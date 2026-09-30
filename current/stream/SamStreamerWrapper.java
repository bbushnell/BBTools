package stream;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.bam.BamIndexWriter;
import structures.ListNum;
import var2.SamFilter;
import var2.ScafMap;
import var2.Scaffold;

/**
 * Command-line wrapper for sequence conversion, SAM filtering and reference-name changes.
 * Uses one consumer with factory-selected readers and writers. Passing records go to
 * the primary output when present; rejected records go to outu when configured, otherwise they
 * are discarded. SAM/BAM input uses raw SamLines unless the primary output requests
 * another sequence format. Non-SAM input uses Reads and disables SAM-specific filters,
 * reference loading, BED filtering, normalization and scaffold renaming.
 * <p>
 * Filters inspect original reference names. When requested, emitted records are renamed
 * before writer handoff, with matching shared SAM headers. Only the raw-SamLine path applies
 * requested CIGAR conversion and normalization, and only to passing records. Rejected
 * records can still be renamed; they are not necessarily byte-for-byte input copies.
 * BAI output delegates directly to BamIndexWriter rather than this streaming pipeline.
 * <p>
 * This is a single-use command object. Parsing and setup change shared I/O, SAM and
 * quality settings without restoring them. Normal streaming waits for writers and
 * checks observed errors, but does not close or join the reader here. Exceptions can
 * bypass normal finalization; this class does not provide a general cleanup guarantee.
 *
 * @author Brian Bushnell, Isla
 * @contributor Shinobu (documentation and formatting)
 * @date November 6, 2025
 */
public class SamStreamerWrapper{

	/*--------------------------------------------------------------*/
	/*----------------            Main              ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs the command and closes a redirected status stream on normal return.
	 * @param args Command-line arguments
	 */
	public static void main(final String[] args){
		final Timer t=new Timer();
		final SamStreamerWrapper x=new SamStreamerWrapper(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses options, checks file access and creates format descriptors without starting
	 * the pipeline Streamer or Writers. Format probing can perform I/O during setup.
	 * Caps shared buffering, enables compression helpers and configures quality/SAM globals.
	 * Disables selected SAM fields for non-SAM or absent primary output unless blocked
	 * by forceparse, CIGAR, BED or secondary-output settings. Retains fields needed by
	 * active contig, position and identity filters. In this optimized branch, CIGAR-dependent
	 * filters disable automatic M-only match derivation; position needs spans, not substitutions.
	 * Identity still resolves M via MD/reference or fails under assertions; CIGARs containing
	 * '=' can still derive match strings. forceparse prevents these changes; it does not restore
	 * globals previously changed elsewhere. Attached-SamLine writing is enabled globally.
	 * @param args Command-line arguments, preparsed for help/configuration and status output
	 */
	SamStreamerWrapper(String[] args){
		{//Preparse block for help, config files, and outstream
			final PreParser pp=new PreParser(args, null/*getClass()*/, false);
			args=pp.args;
			outstream=pp.outstream;
		}

		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();

			threadsIn=parser.threadsIn;
			threadsOut=parser.threadsOut;
			in1=parser.in1;
			out1=parser.out1;
			overwrite=parser.overwrite;
			append=parser.append;
		}

		//Do input/output setup
		fixExtensions();
		checkFileExistence();

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.SAM, null, true, true);
		//STR-041: preserve parsed overwrite/append choices in both preflight and writer descriptors.
		ffout1=FileFormat.testOutput(out1, FileFormat.SAM, null, true, overwrite, append, true);
		ffout2=(out2==null ? null : FileFormat.testOutput(out2, FileFormat.SAM, null, true, overwrite, append, true));

		//Determine if we need to parse SAM fields or can skip for performance.
		//BED filtering needs rname/pos/cigar, and the outu split needs intact records, so both force parsing.
		//STR-042: keep fields required by SamFilter instead of silently changing filter membership.
		if(!forceParse && !fixCigar && !eqx && out2==null && bed==null && (ffout1==null || !ffout1.samOrBam())){
			final boolean filterIdentity=(filter!=null && (filter.minId>0 || filter.maxId<1));
			final boolean filterPosition=(filter!=null && (filter.minPos>Integer.MIN_VALUE || filter.maxPos<Integer.MAX_VALUE));
			final boolean needCigar=(filterIdentity || filterPosition);
			//Identity may resolve M using a reference scaffold; contig membership also needs RNAME.
			SamLine.PARSE_2=(filterIdentity || (filter!=null && filter.contigs!=null));
			SamLine.PARSE_5=needCigar;
			SamLine.PARSE_6=false;
			SamLine.PARSE_7=false;
			SamLine.PARSE_8=false;
			SamLine.PARSE_OPTIONAL=filterIdentity;//Retain MD for ambiguous M identity.
			//SamLine.toRead otherwise demands MD/reference even for position-only M cigars.
			//The '=' exception remains; calcIdentity independently resolves M when identity is requested.
			if(needCigar){SamLine.CONVERT_CIGAR_TO_MATCH=false;}
		}

		//Enable optimizations for SAM->SAM conversion
		ReadStreamByteWriter.USE_ATTACHED_SAMLINE=true;
	}

	/**
	 * Parses wrapper options, SamFilter options, then standard Parser options in argument order.
	 * filter=false disables only SamFilter, not a separately requested BED filter.
	 * Enabling normalization also enables eqx; requesting a SAM version other than 1.3
	 * enables eqx too. Later options can override these flags. Unknown options print a
	 * diagnostic and fail an assertion when assertions are enabled.
	 * @param args Preparsed command-line arguments
	 * @return Parser holding the standard settings copied by the constructor
	 */
	private Parser parse(final String[] args){
		//Create filter
		filter=new SamFilter();
		filter.includeNonPrimary=true;
		filter.includeLengthZero=true;
		filter.includeQfail=true;//Preserve the wrapper's advertised default, independently of SamFilter's default.
		boolean doFilter=true;

		//Create a parser object
		final Parser parser=new Parser();

		//Parse each argument
		for(int i=0; i<args.length; i++){
			final String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			final String[] split=arg.split("=");
			final String a=split[0].toLowerCase();
			final String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("forceparse")){
				forceParse=Parse.parseBoolean(b);
			}else if(a.equals("ref") || a.equals("reference")){
				ref=b;
			}else if(a.equals("rnameasbytes")){
				SamLine.RNAME_AS_BYTES=Parse.parseBoolean(b);
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("samversion") || a.equals("samv") || a.equals("sam")){
				Parser.parseSam(arg, a, b);
				fixCigar=true;
				//SAM 1.4 introduced =/X; other versions enable eqx here. 1.3 leaves the current eqx flag unchanged.
				if(SamLine.VERSION!=1.3f){eqx=true;}
			}else if(a.equals("eqx") || a.equals("eqxcigar") || a.equals("toeqx")){
				//Convert M->=/X (resolving via MD tag or ref=). Synonym for sam=1.4's cigar effect. Default false;
				//INDEPENDENT of normalization (eqx does NOT imply left-shift).
				eqx=Parse.parseBoolean(b);
			}else if(a.equals("normalize") || a.equals("canonicalize") || a.equals("canonicalise")
					|| a.equals("canonize") || a.equals("leftalign")){
				//Deterministic left-alignment of indels vs the reference (alignment canonicalization).
				//Enabling this also enables eqx (later flags can override it); ref= is required, crash-loud if absent.
				normalize=Parse.parseBoolean(b);
				if(normalize){eqx=true;}
			}else if(a.equals("rename") || a.equals("renamescaffolds") || a.equals("renamebylist")){
				//Scaffold rename TSV rewrites SQ header names and record RNAME/RNEXT.
				renameFile=b;
			}else if(a.equals("outu") || a.equals("outunmatched")){
				//Rejected records go here without CIGAR transforms; requested scaffold renaming still applies.
				out2=b;
			}else if(a.equals("bed") || a.equals("bedfile")){
				//BED file for positional filtering (BedReadFilter). Reads are kept/routed by mof+include.
				bed=b;
			}else if(a.equals("minoverlapfraction") || a.equals("mof") || a.equals("overlap")){
				//Fraction of a read's reference span that must lie in the BED to match; 0=any overlap, 1=containment.
				mof=Double.parseDouble(b);
			}else if(a.equals("include") || a.equals("bedinclude") || a.equals("includebed")){
				//true=keep reads matching the BED; false=keep reads that do NOT (reverse/exclude).
				include=Parse.parseBoolean(b);
			}else if(a.equals("filter")){
				doFilter=Parse.parseBoolean(b);
			}else if(filter.parse(arg, a, b)){
				//do nothing

			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(i==0 && !arg.contains("=") && parser.in1==null &&
					FileFormat.isSamOrBamFile(arg) && new File(arg).isFile()){
				parser.in1=arg;
			}else if(i==1 && !arg.contains("=") && parser.out1==null && parser.in1!=null &&
					(FileFormat.isSequence(arg) || FileFormat.isBaiFile(arg))){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		//Disable filter if requested
		if(!doFilter){filter=null;}
		return parser;
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves an existing compressed/uncompressed alternative for the primary input path. */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
	}

	/**
	 * Requires primary input, checks output overwrite/append permissions and input access,
	 * and rejects duplicate primary-input/primary-output/secondary-output paths.
	 * Reference, BED and rename files are loaded later by their respective helpers.
	 */
	private void checkFileExistence(){
		//Ensure input file exists
		if(in1==null){
			throw new RuntimeException("Error - at least one input file is required.");
		}

		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			throw new RuntimeException("\nCan't write to output file "+out1+", "+out2+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read input file "+in1+"\n");
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Primary Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Runs the dedicated BAI branch or constructs and consumes the sequence pipeline once.
	 * Loads requested reference/filter/rename helpers before starting streaming. Shared
	 * headers are retained when either SAM output needs them and, when requested, renamed
	 * before writers start. Reads are requested for non-SAM input or non-SAM primary sequence output;
	 * other streaming routes consume raw SamLines.
	 * Normal completion waits for each writer, captures its counters/error flag, then
	 * checks the reader's observed error flag without calling close or adding a reader
	 * join. Counts printed here are derived from consumed data, not reader aggregates.
	 * Statistics can be printed before the final error exception; they are not proof
	 * of successful completion. Startup/processing exceptions have no common cleanup block.
	 * @param t Timer started before construction; stopped on the normal reporting path
	 */
	void process(final Timer t){
		//Determine processing mode
		final boolean inputSam=(ffin1!=null && ffin1.samOrBam());
		final boolean outputSam=(ffout1!=null && ffout1.samOrBam());
		final boolean outputSam2=(ffout2!=null && ffout2.samOrBam());
		final boolean outputReads=(ffout1!=null && !ffout1.samOrBam());
		final boolean outputBai=(ffout1!=null && ffout1.bai());
		//STR-039: outu needs the input dictionary/metadata even with FASTQ or absent primary output.
		final boolean useSharedHeader=inputSam && (outputSam || outputSam2);
		final boolean makeReads=(outputReads || !inputSam);
		if(!inputSam){
			System.err.println("Input is "+ffin1.formatString()+"; sam filter disabled.");
			filter=null;
			ref=null;
			bed=null;//BED filtering needs reference coordinates, unavailable for non-SAM input
			normalize=false;//no CIGARs to canonicalize for non-SAM input
			renameFile=null;//no scaffold names to rename for non-SAM input
		}

		//Alignment canonicalization needs reference context (rolling can reach beyond the read's footprint)
		if(normalize && ref==null){
			throw new RuntimeException("\nnormalize/canonicalize requires a reference: add ref=<fasta>.\n");
		}

		if(outputBai){
			assert(ffin1.bam()) : "bai output requires bam input.";
			try{
				BamIndexWriter.writeIndex(in1, out1);
			}catch(final Throwable e){
				KillSwitch.exceptionKill(e);
			}
			t.stop();
			outstream.println("Time:                         \t"+t);
			return;
		}

		//Load reference if specified
		if(ref!=null){
			ScafMap.loadReference(ref, true);
			SamLine.RNAME_AS_BYTES=false;
		}

		//Build the BED filter if requested (reads failing it route to outu, or are dropped when outu is unset)
		if(bed!=null){bedFilter=new BedReadFilter(bed, mof, include);}

		//Build the scaffold renamer if requested
		if(renameFile!=null){renamer=new ScaffoldRenamer(renameFile);}

		//Create streamer and writers (fw=primary out; fw2=outu, the non-passing split)
		final Streamer st=StreamerFactory.makeStreamer(ffin1, null, ordered, maxReads, useSharedHeader, makeReads, threadsIn);
		final Writer fw=(ffout1==null ? null : WriterFactory.makeWriter(ffout1, null, threadsOut, null, useSharedHeader));
		final Writer fw2=(ffout2==null ? null : WriterFactory.makeWriter(ffout2, null, threadsOut, null, useSharedHeader));

		//Process data
		st.start();
		//Rename @SQ SN: header lines AFTER the input has loaded the shared header, BEFORE the writer emits it
		//Both raw records and attached SamLines are renamed before writer handoff.
		if(renamer!=null && useSharedHeader){renameSharedHeader();}
		if(fw!=null){fw.start();}
		if(fw2!=null){fw2.start();}

		if(outputReads || !inputSam){
			processAsReads(st, fw, fw2);
		}else{
			processAsSam(st, fw, fw2);
		}

		//Wait for writers to finish
		if(fw!=null){
			errorState|=fw.poisonAndWait();
			readsOut=fw.readsWritten();
			basesOut=fw.basesWritten();
		}
		if(fw2!=null){
			errorState|=fw2.poisonAndWait();
			readsOutu=fw2.readsWritten();
			basesOutu=fw2.basesWritten();
		}

		//Check for errors
		errorState|=st.errorState();

		//Print statistics
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+String.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		outstream.println("Bases Processed:    "+basesProcessed+" \t"+String.format("%.2f Mbp/sec", (basesProcessed/(double)(t.elapsed))*1000));
		if(ffout1!=null){
			outstream.println("Reads Out:          "+readsOut);
			outstream.println("Bases Out:          "+basesOut);
		}
		if(ffout2!=null){
			outstream.println("Reads Outu:         "+readsOutu);
			outstream.println("Bases Outu:         "+basesOutu);
		}

		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException(getClass().getSimpleName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 * Combines the native QC-failure flag with SamFilter before inversion. SamFilter's
	 * includeQfail setting otherwise only affects its external samtools setup (STR-043).
	 * QC failure applies to mapped and unmapped records, before SamFilter's unmapped return.
	 * Disabled filters accept all records; active filters retain SamFilter's null rejection.
	 * @param sl Attached or raw SAM record, possibly null
	 * @return Whether the SAM conditions pass; a separate BED filter may still reject it
	 */
	private boolean passesSamFilter(final SamLine sl){
		if(filter==null){return true;}
		if(sl!=null && !filter.includeQfail && sl.discarded()){return filter.invert;}
		return filter.passesFilter(sl);
	}

	/**
	 * Consumes Read batches until null or an empty input batch, counting data before filtering.
	 * Counts include linked mates; SAM flags alone do not link the converted Read objects.
	 * Filters inspect the primary attachment's original names. Emitted attachments are
	 * renamed once when requested, before either writer receives them; this path performs no explicit
	 * CIGAR conversion or normalization. Containers may be reused, and Read objects are
	 * not copied. Filtered-empty output batches retain input IDs for writer ordering.
	 * @param st Started reader configured to produce Reads
	 * @param fw Primary writer for passing reads, or null
	 * @param fw2 Secondary writer for rejected reads, or null
	 */
	private void processAsReads(final Streamer st, final Writer fw, final Writer fw2){
		final boolean splitting=(fw2!=null);
		for(ListNum<Read> ln=st.nextList(); ln!=null && ln.size()>0; ln=st.nextList()){
			final ArrayList<Read> list=ln.list;
			if(verbose){outstream.println("Got list of size "+ln.size());}

			//Passthrough (out==list, no copy) only when nothing filters or splits the stream
			final boolean passthrough=(filter==null && bedFilter==null && !splitting);
			final ArrayList<Read> out=(passthrough ? list : new ArrayList<Read>(list.size()));
			final ArrayList<Read> outu=(splitting ? new ArrayList<Read>() : null);

			for(final Read r : list){
				//STR-038: non-SAM interleaved input returns R1 with a mate; screen totals count both.
				//Keep reads= in the reader's fragment units while reporting individual reads/bases.
				final int len=r.pairLength();
				readsProcessed+=r.pairCount();
				basesProcessed+=len;

				//A read passes only if it clears every active filter; failures route to outu or are dropped.
				final boolean passes=passesSamFilter(r.samline) &&
					(bedFilter==null || bedFilter.passes(r.samline));
				//STR-040: filter original names, then rename once before either writer receives the batch.
				if(renamer!=null && r.samline!=null && (passes || splitting)){
					renamer.renameRecord(r.samline);
				}
				if(!passes){
					if(outu!=null){outu.add(r);}
					continue;
				}
				if(!passthrough){out.add(r);}
			}

			if(fw!=null){fw.addReads(new ListNum<Read>(out, ln.id));}
			if(fw2!=null){fw2.addReads(new ListNum<Read>(outu, ln.id));}
		}
		if(verbose){outstream.println("Finished.");}
	}

	/**
	 * Consumes raw SAM batches until null or empty, counting each record before filtering.
	 * Passing records receive requested CIGAR operations, optional normalization when
	 * reference bases are available, then optional scaffold renaming. Rejected records skip the
	 * CIGAR operations but are renamed when requested if retained for outu. Filters use original names.
	 * Mutations finish before handoff; original objects and batch IDs are reused, with
	 * new list containers when routing or transformations require them. Empty output
	 * batches still reach their writers. This helper does not finish writers or close input.
	 * @param st Started reader configured to produce SamLines
	 * @param fw Primary writer for passing records, or null
	 * @param fw2 Secondary writer for rejected records, or null
	 */
	private void processAsSam(final Streamer st, final Writer fw, final Writer fw2){
		final boolean splitting=(fw2!=null);
		final boolean transforming=(fixCigar || eqx);
		final ScafMap scafMap=(normalize ? ScafMap.defaultScafMap() : null);//reference for indel left-alignment
		for(ListNum<SamLine> ln=st.nextLines(); ln!=null && ln.size()>0; ln=st.nextLines()){
			final ArrayList<SamLine> list=ln.list;
			if(verbose){outstream.println("Got list of size "+ln.size());}

			//Passthrough (out==list, no copy) only when nothing filters, transforms, renames, or splits the stream
			final boolean passthrough=(filter==null && bedFilter==null && renamer==null && !transforming && !splitting);
			final ArrayList<SamLine> out=(passthrough ? list : new ArrayList<SamLine>(list.size()));
			final ArrayList<SamLine> outu=(splitting ? new ArrayList<SamLine>() : null);

			for(final SamLine sl : list){
				final int len=sl.lengthOrZero();
				readsProcessed++;
				basesProcessed+=len;

				//A read passes only if it clears every active filter (SamFilter AND BedReadFilter).
				//Rejected records skip CIGAR transforms but can be renamed for outu, or dropped if outu is unset.
				final boolean passes=passesSamFilter(sl) &&
					(bedFilter==null || bedFilter.passes(sl));//filters use ORIGINAL scaffold names
				if(!passes){
					if(outu!=null){
						if(renamer!=null){renamer.renameRecord(sl);} //outu records must match the renamed header too
						outu.add(sl);
					}
					continue;
				}

				//Length-zero guard for non-SAM output (kept reads only); unchanged behavior for SAM/BAM output
				if(!(len>0 || ffout1==null || ffout1.samOrBam())){continue;}

				//CIGAR transforms on KEPT reads. eqx (convert M->=/X) and sam-version conversion are INDEPENDENT.
				if(transforming && sl.cigar!=null){
					if(eqx){
						//eqx (= sam=1.4): emit =/X, resolving M via the MD tag or loaded reference (ScafMap from ref=).
						//Crash-loud under -ea on an M-only cigar with neither MD nor ref: do what was asked, or fail -- never
						//silently keep an unconverted cigar or drop it.
						rebuildCigar14(sl, false);
					}else{
						sl.setCigar(SamLine.toCigar13(sl.cigar));//Version conversion with eqx disabled: =/X->M
					}
				}

				//Left-align indels vs the reference after any requested eqx conversion.
				if(normalize && sl.cigar!=null){
					final Scaffold scaf=scafMap.getScaffold(sl.rnameS());
					if(scaf!=null && scaf.bases!=null){CigarNormalizer.normalize(sl, scaf.bases);}
				}

				//Scaffold rename LAST, after BED filter + cigar transforms (which all use the original names)
				if(renamer!=null){renamer.renameRecord(sl);}

				if(!passthrough){out.add(sl);}
			}

			if(fw!=null){fw.addLines(new ListNum<SamLine>(out, ln.id));}
			if(fw2!=null){fw2.addLines(new ListNum<SamLine>(outu, ln.id));}
		}
		if(verbose){outstream.println("Finished.");}
	}

	/** Requests the shared header with waiting enabled and replaces renamed SQ entries in its mutable list.
	 * Called after input startup but before either writer starts; both SAM outputs use
	 * the resulting names. A null header is left unchanged; no defensive list copy is made.
	 */
	private void renameSharedHeader(){
		final ArrayList<byte[]> hdr=SamReadInputStream.getSharedHeader(true);
		if(hdr==null){return;}
		for(int i=0; i<hdr.size(); i++){
			final String line=new String(hdr.get(i));
			final String renamed=renamer.renameHeaderLine(line);
			if(!renamed.equals(line)){hdr.set(i, renamed.getBytes());}
		}
	}

	/** Rebuilds a SamLine's CIGAR in SAM 1.4 (=/X) representation from its match string.
	 * @param sl Non-null record whose CIGAR will be replaced
	 * @param allowM true keeps M ops (version-conversion only); false resolves M->=/X via MD tag or
	 * loaded reference; with assertions enabled, unresolved M without either causes a loud failure.
	 */
	private static void rebuildCigar14(final SamLine sl, final boolean allowM){
		final byte[] shortMatch=sl.toShortMatch(allowM);
		final byte[] longMatch=Read.toLongMatchString(shortMatch);
		final int start=sl.pos-1;
		final int stop=start+Read.calcMatchLength(longMatch)-1;
		sl.setCigar(SamLine.toCigar14(longMatch, start, stop, Integer.MAX_VALUE, sl.seq));
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** SAM/BAM filter for quality/mapping criteria; null disables filtering. */
	private SamFilter filter;
	/** Positional BED-overlap filter; null disables BED filtering. Built in process() from bed/mof/include. */
	private BedReadFilter bedFilter;
	/** Scaffold renamer (old->new); null disables renaming. Built in process() from renameFile. */
	private ScaffoldRenamer renamer;

	/** Single primary input path; SAM is the fallback format. */
	private String in1=null;
	/** Primary sequence output path, BAI output path, or null for no primary output. */
	private String out1=null;
	/** Secondary output for rejected records, which skip CIGAR transforms but can be renamed; null disables it. */
	private String out2=null;
	/** Permission to replace existing output files when not appending. */
	private boolean overwrite=true;
	/** Append flag for sequence writers; also permits existing output paths when overwrite is false. */
	private boolean append=false;
	/** Reference path for coordinate lookups or CIGAR reconstruction. */
	private String ref=null;
	/** BED file for positional filtering; null disables it. */
	private String bed=null;
	/** Two-column old/new scaffold-name TSV path; null disables renaming. */
	private String renameFile=null;
	/** Minimum overlap fraction for the BED filter (0=any overlap, 1=containment). */
	private double mof=0.0;
	/** BED filter sense: true keeps reads matching the BED, false keeps non-matching. */
	private boolean include=true;

	/** Input file format descriptor. */
	private FileFormat ffin1;
	/** Output file format descriptor. */
	private FileFormat ffout1;
	/** Secondary (outu) output file format descriptor; null when no split is requested. */
	private FileFormat ffout2;

	/** Reader-selection thread hint; negative selects the format default. */
	private int threadsIn=-1;
	/** Writer-selection thread hint; negative selects the format default. */
	private int threadsOut=-1;

	/** Individual reads counted before filtering, including linked mates in Read mode. */
	private long readsProcessed=0;
	/** Primary writer's reported read count after its normal completion. */
	private long readsOut=0;
	/** Secondary writer's reported read count after its normal completion. */
	private long readsOutu=0;
	/** Bases counted before filtering, including linked mate bases in Read mode. */
	private long basesProcessed=0;
	/** Primary writer's reported base count after its normal completion. */
	private long basesOut=0;
	/** Secondary writer's reported base count after its normal completion. */
	private long basesOutu=0;

	/*--------------------------------------------------------------*/

	/** Accumulates observed writer and reader errors on the normal finalization path. */
	public boolean errorState=false;
	/** Ordering request passed to the input factory; output descriptors request ordering. */
	public boolean ordered=true;
	/** Input limit forwarded to the selected reader in its units, before filters; negative is unlimited. */
	private long maxReads=-1;
	/**
	 * Prevents this constructor from disabling SAM parse flags; does not restore prior global settings.
	 */
	private boolean forceParse;
	/** Enables raw-route CIGAR version conversion; eqx selects whether to resolve M or collapse =/X. */
	private boolean fixCigar;
	/** Resolve raw-route M CIGAR operations via MD/reference; also set by version/normalization options. */
	private boolean eqx=false;
	/** Requests raw-route indel normalization with reference bases; enabling it also sets eqx during parsing. */
	private boolean normalize=false;

	/*--------------------------------------------------------------*/

	/** Output stream used for status messages and timing summaries. */
	private PrintStream outstream=System.err;
	/** Enables verbose debugging output. */
	public static boolean verbose=false;

}
