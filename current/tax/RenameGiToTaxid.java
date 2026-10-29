package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.HashArray1D;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import structures.IntList;
import tracker.ReadStats;

/**
 * Renames FASTA or GFF sequence identifiers to taxonomy IDs.
 * <p>
 * Local mode resolves complete headers against a loaded {@link TaxTree}; server
 * mode batches accessions, GI numbers, UNITE identifiers or full headers through
 * {@link TaxClient} and rewrites FASTA headers, GFF seqids and matching
 * {@code ##sequence-region} pragmas. GFF output always preserves the original
 * feature columns, so GFF input rejects {@code prefix=f}.
 *
 * @author Brian Bushnell
 * @date Mar 10, 2015
 */
public class RenameGiToTaxid{

	/** Runs taxonomy-ID header renaming with the supplied command-line options.
	 * @param args BBTools-style {@code flag=value} arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		RenameGiToTaxid x=new RenameGiToTaxid(args);
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/** Parses arguments, configures I/O formats, and loads local taxonomy tables or the remote-server route.
	 * @param args BBTools-style {@code flag=value} arguments */
	public RenameGiToTaxid(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_BGZIP=ReadWrite.USE_UNBGZIP=ReadWrite.PREFER_BGZIP=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("prefix")){
				prefix=Parse.parseBoolean(b);

			}else if(a.equals("server") || a.equals("useserver")){
				if(b!=null && b.startsWith("http")){
					useServer=true;
					String path=b;
					if(!path.endsWith("/")){path+="/";}
					Shared.setTaxServer(path);
				}else{
					useServer=Parse.parseBoolean(b);
				}
			}else if(a.equals("title")){
				title=(b==null ? ">" : (">"+b+"|")).getBytes();
				//TODO: Probable bug - invalidTitle is not updated to match a custom title.
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
			}else if(a.equals("accession")){
				accessionFile=b;
			}else if(a.equals("pattern")){
				patternFile=b;
			}else if(a.equals("tree") || a.equals("taxtree")){
				useTree=TaxTree.parseTreeFlag(b);
			}else if(a.equals("invalid")){
				outInvalid=b;
			}else if(a.equals("deleteinvalid")){
				deleteInvalid=Parse.parseBoolean(b);
			}else if(a.equals("badheaders")){
				badHeaders=b;
			}else if(a.equals("maxbadheaders") || a.equals("maxinvalidheaders")){
				maxInvalidHeaders=Parse.parseKMG(b);
			}else if(a.equals("keepall")){
				keepAll=Parse.parseBoolean(b);
			}else if(a.equals("shrinknames")){
				shrinkNames=Parse.parseBoolean(b);
			}else if(a.equals("warn")){
				warnBadHeaders=Parse.parseBoolean(b);
			}

			else if(a.equals("maxpigzprocesses")){
				AccessionToTaxid.maxPigzProcesses=Integer.parseInt(b);
			}else if(a.equals("skipparse")){
				AccessionToTaxid.skipParse=Parse.parseBoolean(b);
			}else if(a.equals("skiphash")){
				AccessionToTaxid.skipHash=Parse.parseBoolean(b);
			}

			else if(a.equals("mode")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					mode=Integer.parseInt(b);
				}else if("accession".equalsIgnoreCase(b)){
					mode=ACCESSION_MODE;
				}else if("unite".equalsIgnoreCase(b)){
					mode=UNITE_MODE;
					TaxTree.UNITE_MODE=true;
				}else if("gi".equalsIgnoreCase(b)){
					mode=GI_MODE;
				}else if("header".equalsIgnoreCase(b)){
					mode=HEADER_MODE;
				}else{
					assert(false) : "Bad mode: "+b;
				}
			}

			else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("in") || a.equals("in1")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(new File(b).exists()){
					in1.add(b);
				}else{
					for(String bb : b.split(",")){
						in1.add(bb);
					}
				}
			}else if(new File(arg).exists()){ //For asterisk expansion
				in1.add(arg);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		if(useServer){
			giTableFile=null;
			accessionFile=null;
			patternFile=null;
			if(mode!=UNITE_MODE){useTree=false;}
		}//else if taxpath!=null... set them

		{//Process parser fields
			Parser.processQuality();

			maxReads=parser.maxReads;

			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out1=parser.out1;
		}

		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}
		if("auto".equalsIgnoreCase(accessionFile)){accessionFile=TaxTree.defaultAccessionFile();}
		if("auto".equalsIgnoreCase(patternFile)){patternFile=TaxTree.defaultPatternFile();}

		assert(FastaReadInputStream.settingsOK());

		if(in1==null || in1.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		assert(out1!=null) : "This program requires an output file.";

		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");
		}

		ffout1=FileFormat.testOutput(out1, FileFormat.FA, null, true, overwrite, append, false);
		ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.FA, null, true, overwrite, append, false);
		ffin1=new ArrayList<FileFormat>(in1.size());
		for(String s : in1){
			FileFormat ff=FileFormat.testInput(s, FileFormat.FA, null, true, true);
			ffin1.add(ff);
		}
		if(!prefix){
			for(FileFormat ff : ffin1){
				if(ff.gff()){throw new RuntimeException(GFF_PREFIX_FALSE_MESSAGE);}
			}
		}

		if(ffoutInvalid!=null){keepAll=false;}

		assert(giTableFile!=null || accessionFile!=null || TaxTree.SILVA_MODE || useServer) : "No gi or accession information loaded.";

		if(useTree){
			tree=TaxTree.loadTaxTree(outstream, true, false);
			assert(tree.nameMap!=null);
		}else{
			tree=null;
			if(!useServer){throw new RuntimeException("No tree specified.");}
		}

		if(giTableFile!=null){
			GiToTaxid.initialize(giTableFile);
		}

		if(patternFile!=null){
			Timer t=new Timer();
			AnalyzeAccession.loadCodeMap(patternFile);
			outstream.println("Loading pattern table.");
			t.stopAndPrint();
		}

		if(accessionFile!=null){
			AccessionToTaxid.tree=tree;
			outstream.println("Loading accession table.");
			AccessionToTaxid.load(accessionFile);
		}
	}

	/** Renames every input into the primary and optional invalid-header outputs, then reports FASTA or GFF totals.
	 * @param t active timer for the processing phase */
	void process(Timer t){

		ByteStreamWriter bsw=(ffout1==null ? null : new ByteStreamWriter(ffout1)); //Actually, this is required.
		if(bsw!=null){bsw.start();}

		ByteStreamWriter bswInvalid=null;
		if(ffoutInvalid!=null){
			bswInvalid=new ByteStreamWriter(ffoutInvalid);
			bswInvalid.start();
		}

		ByteStreamWriter bswBadHeaders=null;
		if(badHeaders!=null){
			bswBadHeaders=new ByteStreamWriter(badHeaders, overwrite, append, false);
			bswBadHeaders.start();
		}

		final HashArray1D counts=(countTable && !prefix) ? new HashArray1D(256000, -1L, true) : null;

		gffIn=false;
		for(FileFormat ffin : ffin1){
			//gffIn remains set after the first GFF input, selecting final line/byte reporting
			//and GFF classification for later server-mode batches globally.
			gffIn=gffIn||ffin.gff();
			ByteFile bf=ByteFile.makeByteFile(ffin);
			if(useServer){
				processInner_server(bf, bsw, bswInvalid, bswBadHeaders, counts, ffin.format());
			}else{
				processInner(bf, bsw, bswInvalid, bswBadHeaders, counts, null);
			}
		}

		if(bsw!=null){
			errorState|=bsw.poisonAndWait();
			if(deleteInvalid && invalidReads>0 && !ffout1.stdio()){
				try{
					System.err.println("Deleting "+out1);
					new File(out1).delete();
				}catch(Exception e){
					System.err.println("An error occured while attempting to delete "+out1);
					e.printStackTrace();
				}
			}
		}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}
		if(bswBadHeaders!=null){errorState|=bswBadHeaders.poisonAndWait();}

		t.stop();
		if(!gffIn){
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));

			outstream.println();
			outstream.println("Valid Sequences:   \t"+validReads);
			outstream.println("Valid Bases:       \t"+validBases);
			outstream.println("Invalid Sequences: \t"+invalidReads);
			outstream.println("Invalid Bases:     \t"+invalidBases);
		}else{
			outstream.println(Tools.timeLinesBytesProcessed(t, linesIn, basesProcessed, 8));

			outstream.println();
			outstream.println("Valid Lines:       \t"+validLines);
			outstream.println("Valid Bytes:       \t"+validBases);
			outstream.println("Invalid Lines:     \t"+invalidLines);
			outstream.println("Invalid Bytes:     \t"+invalidBases);
		}
		if(counts!=null){
			outstream.println("Unique Taxa:       \t"+taxaCounted);
		}

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}


	/** Processes one local FASTA stream, resolving each header against the loaded {@link TaxTree}.
	 * @param bf input file reader, closed on normal return
	 * @param bsw primary output stream receiving valid records and, when {@link #keepAll} is true, invalid records
	 * @param bswInvalid optional output stream receiving invalid records when {@link #keepAll} is false
	 * @param bswBadHeaders optional writer for raw headers that cannot be resolved
	 * @param counts optional per-taxID occurrence counter used only when prefix output is disabled
	 * @param ids optional precomputed taxID list aligned to FASTA headers */
	private void processInner(ByteFile bf, ByteStreamWriter bsw, ByteStreamWriter bswInvalid, ByteStreamWriter bswBadHeaders, HashArray1D counts, IntList ids){

		int readsProcessedInner=0;

		byte[] line=bf.nextLine();
		boolean valid=false;
		while(line!=null){
			if(line.length>0 && line[0]=='>'){
				readsProcessedInner++;
				readsProcessed++;
				if(maxReads>0 && readsProcessed>maxReads){break;}
				int initial=1, terminal=line.length;
				final int number;
				if(ids==null){
					final TaxNode tn;

					{
						{
							//					Handles renumbering when the format is correct but the number is wrong.
							if(Tools.startsWith(line, ">tid|")){
								initial=6;
								while(initial<=line.length && line[initial-1]!='|'){initial++;}
							}else if(Tools.startsWith(line, ">ncbi|")){
								initial=7;
								while(initial<=line.length && line[initial-1]!='|'){initial++;}
							}
						}
						// FIXED [tax/RenameGiToTaxid#001]: for ">tid|123" or
						// ">ncbi|123" with no second '|', the seek loop advanced
						// initial past the end and then sliced a negative-length
						// header suffix. Decouple lookupStart from the preserved
						// suffix start: preserve nothing, then look up the whole
						// tid/ncbi token so parseNodeFromHeader resolves it.
						final int lookupStart;
						if(initial>line.length){initial=terminal; lookupStart=1;}
						else{lookupStart=initial;}

						if(shrinkNames){//This is for nr/nt
							for(int i=initial; i<terminal; i++){
								if(line[i]==1){//SOH
									terminal=i;
								}
							}
						}

						String s=new String(line, lookupStart, terminal-lookupStart);

						tn=tree.parseNodeFromHeader(s, true);
					}
					number=(tn==null ? -1 : tn.id);
				}else{
					number=ids.get((int)(readsProcessedInner-1));

					if(shrinkNames){//This is for nr/nt
						for(int i=initial; i<terminal; i++){
							if(line[i]==1){//SOH
								terminal=i;
							}
						}
					}
				}

				valid=(number>=0);
				if(valid){
					validReads++;
					bsw.print(title);
					bsw.print(number);
					if(prefix){
						bsw.print('|');
						for(int i=initial; i<terminal; i++){
							bsw.print(line[i]);
						}
					}else if(counts!=null){
						bsw.print('|');
						int count=counts.increment(number, 1);
						bsw.print(count);
						if(count==1){taxaCounted++;}
					}
					bsw.println();
				}else{
					invalidReads++;
					if(deleteInvalid){
						System.err.println("Invalid sequence detected; aborting.\n");
						break;
					}
					if(bswBadHeaders!=null){bswBadHeaders.println(line);}
					if(maxInvalidHeaders>=0 && invalidReads>maxInvalidHeaders){
						KillSwitch.kill("Maximum bad headers exceeded: "+maxInvalidHeaders+"\n"+new String(line));
					}
					if(keepAll){
						if(shrinkNames){
							for(int i=0; i<terminal; i++){
								bsw.print(line[i]);
							}
							bsw.println();
						}else{
							bsw.println(line);
						}
					}else if(bswInvalid!=null){
						if(shrinkNames){
							for(int i=0; i<terminal; i++){
								bswInvalid.print(line[i]);
							}
							bswInvalid.println();
						}else{
							bswInvalid.println(line);
						}
					}
				}
			}else{
				basesProcessed+=line.length;
				if(valid || keepAll){
					if(valid){validBases+=line.length;}
					else{invalidBases+=line.length;}
					bsw.println(line);
				}else{
					invalidBases+=line.length;
					if(bswInvalid!=null){
						bswInvalid.println(line);
					}
				}
			}
			line=bf.nextLine();
		}

		errorState|=bf.close();
	}

	/** Returns whether a bad header's leading token looks like a plausible accession for warning output.
	 * @param line header or GFF feature line to evaluate
	 * @return {@code true} if the first token is 4-18 accession-like characters with an optional version suffix */
	private static boolean looksLikeRealAccession(byte[] line){
		int space=Tools.indexOf(line, ' ');
		if(space<0){space=line.length;}
		if(space>18 || space<4){return false;}
		//...  hmm...  this is a pretty short list for false cases!
		int dot=-1;
		for(int i=0; i<space; i++){
			if(line[i]=='.'){
				if(dot>=0){return false;}//Only 1 dot allowed
				dot=i;
			}
		}
		if(dot>0){
			if(dot!=space-2){return false;}
		}
		for(int i=0; i<space; i++){
			byte b=line[i];
			if(b!='_' && b!='-' && b!='.' && !Tools.isLetterOrDigit(b)){return false;}
		}
		return true;
	}

	/** Appends the mode-selected server query token from a FASTA/FASTQ header, followed by a comma.
	 * @param line header line beginning with {@code >} or {@code @}
	 * @param bb destination for the extracted identifier and delimiter */
	void appendHeaderLine(byte[] line, ByteBuilder bb){
		assert(line[0]=='>' || line[0]=='@') : new String(line);

		if(mode==ACCESSION_MODE){
			for(int i=1; i<line.length; i++){
				byte b=line[i];
				if(b==' ' || b=='.'){break;}
				else{bb.append(b);}
			}
		}else if(mode==GI_MODE){
			for(int i=1; i<line.length; i++){
				byte b=line[i];
				if(b==' ' || b=='|'){break;}
				else{bb.append(b);}
			}
		}else if(mode==UNITE_MODE){
			int initial=Tools.indexOf(line, '|');
			for(int i=initial+1; i<line.length; i++){
				byte b=line[i];
				if(b==' ' || b=='.' || b=='|'){break;}
				else{bb.append(b);}
			}
		}else{
			for(int i=1; i<line.length; i++){
				byte b=line[i];
				bb.append(b);
			}
		}
		bb.append(',');
	}

	/** Rewrites one buffered server-mode batch in place using the handler for its input format.
	 * @param lines buffered FASTA or GFF lines
	 * @param counts optional per-taxID occurrence counter used only when prefix output is disabled
	 * @param bswBadHeaders optional writer for identifiers that cannot be resolved
	 * @param format {@link FileFormat#FA} or {@link FileFormat#GFF} */
	private void updateHeadersFromServer(ArrayList<byte[]> lines, HashArray1D counts, ByteStreamWriter bswBadHeaders, int format){
		if(format==FileFormat.FA){
			updateHeadersFromServer_fasta(lines, counts, bswBadHeaders);
		}else if(format==FileFormat.GFF){
			updateHeadersFromServer_gff(lines, counts, bswBadHeaders);
		}else{
			assert(false) : "Unsupported type: "+format;
		}
	}

	/** Replaces non-{@code >tid} FASTA headers in one server-mode batch.
	 * @param lines FASTA lines to rewrite in place
	 * @param counts optional per-taxID occurrence counter used only when prefix output is disabled
	 * @param bswBadHeaders optional writer for unresolvable headers */
	private void updateHeadersFromServer_fasta(ArrayList<byte[]> lines, HashArray1D counts, ByteStreamWriter bswBadHeaders){
		ByteBuilder bb=new ByteBuilder();
		ArrayList<String> names=new ArrayList<String>();
		int terms=0;
		for(byte[] line : lines){
			if(line[0]=='>' && !Tools.startsWith(line, ">tid")){
				terms++;
				appendHeaderLine(line, bb);//At this point bb contains something like "GHE08303,ABC123,KR99,"
				if(mode==UNITE_MODE){
					int bar=Tools.indexOf(line, '|');
					names.add(new String(line, 1, bar-1));
				}
			}
		}
		if(bb.length()<1){return;}

		assert(terms>0) : bb.count(',')+", "+terms+"\n"+bb+"\n";
		assert(bb.endsWith(',')) : bb.count(',')+", "+terms+"\n"+bb+"\n";
		assert(bb.count(',')==terms) : bb.count(',')+", "+terms+"\n"+bb+"\n";
		bb.length--;//At this point bb contains something like "GHE08303,ABC123,KR99"


		final int[] serverIds=translateToTaxIDs(bb, terms);
		assert(serverIds!=null) : "Null response for '"+bb.toString()+"'";
		assert(serverIds.length==terms) : "Incorrect response length for query."
				+"\nQuery size "+terms+":\n"+bb
				+"\nResponse size "+serverIds.length+":\n"+Arrays.toString(serverIds)+"\n";
		bb.clear();

		if(!names.isEmpty()){
			assert(tree!=null) : "Need to load a TaxTree.";
			assert(names.size()==serverIds.length);
			for(int i=0; i<serverIds.length; i++){
				final String name=names.get(i);
				if(serverIds[i]<0){
					TaxNode tn=tree.getNodeByName(name);
					if(tn!=null){serverIds[i]=tn.id;}
				}else{
					//Sometimes the species gets renamed.
				}
			}
		}

		for(int lineNum=0, serverNum=0; lineNum<=lines.size(); lineNum++){
			byte[] line=lines.get(lineNum);
			if(line[0]=='>' && !Tools.startsWith(line, ">tid")){
				bb.clear();
				final int tid=serverIds[serverNum];
				if(tid<0){
					//WARN
					if(bswBadHeaders!=null){
						bswBadHeaders.print(tid).tab();
						bswBadHeaders.print(looksLikeRealAccession(line)).tab();
						bswBadHeaders.println(line);
					}else if(warnBadHeaders){
						System.err.println(tid+"\t"+looksLikeRealAccession(line)+"\t"+new String(line));
					}
				}
				int initial=1, terminal=line.length;
				if(shrinkNames){//This is for nr/nt
					for(int i=initial; i<terminal; i++){
						if(line[i]==1){//SOH
							terminal=i;
						}
					}
				}

				bb.append(title);
				bb.append(tid);
				if(prefix){
					bb.append('|');
					for(int i=initial; i<terminal; i++){
						bb.append(line[i]);
					}
				}else if(counts!=null && tid>=0){
					bb.append('|');
					int count=counts.increment(tid, 1);
					bb.append(count);
					if(count==1){taxaCounted++;}
				}

				lines.set(lineNum, bb.toBytes());

				serverNum++;
				if(serverNum>=serverIds.length){break;}
			}
		}
		if(maxInvalidHeaders>=0 && invalidReads>maxInvalidHeaders){
			//TODO: Probable bug - the current batch's invalid reads are counted later in dumpBuffer.
			KillSwitch.kill("Maximum bad headers exceeded: "+maxInvalidHeaders);
		}
	}

	/** Translates comma-separated identifiers through the selected server endpoint, retrying null or wrong-length replies.
	 * @param bb comma-separated query terms
	 * @param terms expected number of taxonomy IDs in the reply
	 * @return aligned taxonomy IDs, or null after all retries fail */
	private int[] translateToTaxIDs(ByteBuilder bb, int terms){
		//This block was created due to mysterious crashes processing nr.faa.
		//They seemed to be due to incomplete responses.
		int[] response=null;
		final String query=bb.toString();
		for(int i=0; i<10 && response==null; i++){
			if(mode==ACCESSION_MODE || mode==UNITE_MODE){
				response=TaxClient.accessionToTaxidArray(query);
			}else if(mode==GI_MODE){
				response=TaxClient.giToTaxidArray(query);
			}else{
				response=TaxClient.headerToTaxidArray(query);//TODO: Keep this as a ByteBuilder
			}
			if(response==null){
				System.err.println("Retrying null response for length "+terms+" query: "+query);
			}else if(response.length!=terms){
				System.err.println("Retrying corrupt length "+response.length+" response for length "+terms+" query: "+query);
				response=null;
			}
			if(response==null){
				try{
					Thread.sleep(500*(i+2));
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
		}
		return response;
	}

	/** Replaces GFF feature seqids and matching {@code ##sequence-region} seqids in one server-mode batch.
	 * @param lines GFF lines to rewrite in place
	 * @param counts optional per-taxID occurrence counter used only when prefix output is disabled
	 * @param bswBadHeaders optional writer for unresolvable GFF seqids */
	private void updateHeadersFromServer_gff(ArrayList<byte[]> lines, HashArray1D counts, ByteStreamWriter bswBadHeaders){
		assert(prefix) : "Constructor rejects GFF input with prefix=f.";
		ByteBuilder bb=new ByteBuilder();
		int terms=0;
		//[tax/RenameGiToTaxid#004] Query seqids from BOTH feature lines (column 1) AND ##sequence-region
		//pragmas, so the pragma seqids get renamed identically to the features (and thus to the paired fna
		//headers). Without this, a renamed gff keeps stale accession-named ##sequence-region lines that
		//disagree with its own features and with the fna, breaking any tool that reads the pragma's seqid.
		for(byte[] line : lines){
			if(line.length==0){continue;}
			if(Tools.startsWith(line, SEQREGION_PRAGMA)){
				String sid=sequenceRegionId(line);
				if(sid!=null && !sid.startsWith("tid|")){
					terms++;
					if(bb.length()>0){bb.append(',');}
					bb.append(sid);
				}
			}else if(line[0]!='#' && !Tools.startsWith(line, "tid")){
				terms++;
				if(bb.length()>0){bb.append(',');}
				for(byte b : line){
					if(b=='\t'){break;}
					bb.append(b);
				}
			}
		}
		if(bb.length()<1){return;}

		// HARDENED [tax/RenameGiToTaxid#002]: a wrong-length GFF server
		// response used to silently misalign subsequent seqids because only
		// null was retried. FASTA already routes through translateToTaxIDs;
		// use the same retry path for GFF, with persistent failure still
		// returning null so the KillSwitch fires loudly.
		int[] serverIds=translateToTaxIDs(bb, terms);
		if(serverIds==null){
			KillSwitch.kill("Null response for '"+bb.toString()+"'");
		}
		bb.clear();

		//Rewrite feature lines and ##sequence-region pragmas in the SAME order they were queried, so
		//serverNum stays aligned to serverIds. Loop bound is strict (< size) and gated on serverNum so it
		//can never index past either array.
		for(int lineNum=0, serverNum=0; lineNum<lines.size() && serverNum<serverIds.length; lineNum++){
			byte[] line=lines.get(lineNum);
			if(line.length==0){continue;}
			if(Tools.startsWith(line, SEQREGION_PRAGMA)){
				String sid=sequenceRegionId(line);
				if(sid!=null && !sid.startsWith("tid|")){
					lines.set(lineNum, renameSeqRegion(line, sid, serverIds[serverNum]));
					serverNum++;
				}
			}else if(line[0]!='#' && !Tools.startsWith(line, "tid")){
				bb.clear();
				final int tid=serverIds[serverNum];
				if(tid<0){
					//WARN
					if(bswBadHeaders!=null){
						bswBadHeaders.print(tid).tab();
						bswBadHeaders.print(looksLikeRealAccession(line)).tab();
						bswBadHeaders.println(line);
					}else if(warnBadHeaders){
						System.err.println(tid+"\t"+looksLikeRealAccession(line)+"\t"+new String(line));
					}
				}

				bb.append("tid|");
				bb.append(tid);
				if(prefix){
					bb.append('|');
					bb.append(line);
				}else if(counts!=null && tid>=0){
					bb.append('|');
					int count=counts.increment(tid, 1);
					bb.append(count);
					if(count==1){taxaCounted++;}
				}

				lines.set(lineNum, bb.toBytes());
				serverNum++;
			}
		}
		if(maxInvalidHeaders>=0 && invalidReads>maxInvalidHeaders){
			//TODO: Probable bug - the current batch's invalid features are counted later in dumpBuffer.
			KillSwitch.kill("Maximum bad headers exceeded: "+maxInvalidHeaders);
		}
	}

	/** Returns the seqid from a {@code ##sequence-region seqid start end} pragma line, or null if absent.
	 * @param line raw GFF pragma line
	 * @return parsed seqid, or null if the line has no seqid */
	private static String sequenceRegionId(byte[] line){
		int i=SEQREGION_PRAGMA.length();
		while(i<line.length && (line[i]==' ' || line[i]=='\t')){i++;}
		final int a=i;
		while(i<line.length && line[i]!=' ' && line[i]!='\t'){i++;}
		return i>a ? new String(line, a, i-a) : null;
	}

	/** Rewrites a {@code ##sequence-region} line's seqid to {@code tid|tid|seqid}
	 * (or {@code tid|tid} when {@link #prefix} is false), preserving the pragma
	 * prefix and trailing start/end coordinates.
	 * @param line raw GFF pragma line
	 * @param sid seqid parsed from the pragma
	 * @param tid taxonomy ID assigned to the seqid
	 * @return renamed GFF pragma line */
	private byte[] renameSeqRegion(byte[] line, String sid, int tid){
		int i=SEQREGION_PRAGMA.length();
		while(i<line.length && (line[i]==' ' || line[i]=='\t')){i++;}
		final int a=i;                       //seqid start
		final int b=a+sid.length();          //seqid end (first whitespace after it)
		ByteBuilder bb=new ByteBuilder(line.length+16);
		bb.append(line, 0, a);               //"##sequence-region" + its trailing spacing
		bb.append("tid|").append(tid);
		if(prefix){bb.append('|').append(sid);}
		bb.append(line, b, line.length-b);   //trailing start/end coordinates
		return bb.toBytes();
	}

	/** Processes one server-mode FASTA or GFF stream in triggered byte batches.
	 * @param bf input file reader, closed on normal return
	 * @param bsw primary output stream receiving valid records and, when {@link #keepAll} is true, invalid records
	 * @param bswInvalid optional output stream receiving invalid records when {@link #keepAll} is false
	 * @param bswBadHeaders optional writer for identifiers that cannot be resolved
	 * @param counts optional per-taxID occurrence counter used only when prefix output is disabled
	 * @param format {@link FileFormat#FA} or {@link FileFormat#GFF} */
	private void processInner_server(ByteFile bf, ByteStreamWriter bsw, ByteStreamWriter bswInvalid, ByteStreamWriter bswBadHeaders, HashArray1D counts, int format){

		ArrayList<byte[]> lines=new ArrayList<byte[]>();
		byte[] line=bf.nextLine();
		boolean valid=false;
		long storedBytes=0;
		long readsBuffered=readsProcessed;

		while(line!=null){

			if(line.length>0){
				if(format==FileFormat.FA && line[0]=='>'){
					if(maxReads>0 && readsBuffered>=maxReads){break;}
					readsBuffered++;
				}
				linesIn++;
				lines.add(line);
				storedBytes+=line.length;
				if(storedBytes>=maxStoredBytes){
					updateHeadersFromServer(lines, counts, bswBadHeaders, format);
					valid=dumpBuffer(lines, valid, bsw, bswInvalid);
					lines=new ArrayList<byte[]>();
					storedBytes=0;
					if(deleteInvalid && invalidReads>0){
						System.err.println("Invalid sequence detected; aborting.\n"
								+ "Input file:  \t"+bf.name()+"\n"
								+ "Output file: \t"+(bsw==null ? "null" : bsw.fname)+"\n"
								+ "Line:        \t"+new String(line)+"\n");
						break;
					}
				}
			}
			line=bf.nextLine();
		}

		if(storedBytes>0){
			updateHeadersFromServer(lines, counts, bswBadHeaders, format);
			valid=dumpBuffer(lines, valid, bsw, bswInvalid);
			lines=new ArrayList<byte[]>();
			storedBytes=0;
		}

		errorState|=bf.close();
	}

	/** Emits a rewritten server-mode batch and updates valid/invalid FASTA or GFF counters.
	 * @param lines rewritten batch to emit
	 * @param valid validity state carried from the previous FASTA batch
	 * @param bsw primary output stream receiving valid records and, when {@link #keepAll} is true, invalid records
	 * @param bswInvalid optional output stream receiving invalid records when {@link #keepAll} is false
	 * @return validity state to carry to the next FASTA batch */
	private boolean dumpBuffer(ArrayList<byte[]> lines, boolean valid, ByteStreamWriter bsw, ByteStreamWriter bswInvalid){

		for(byte[] line : lines){

			if(line.length>0 && line[0]=='>'){
				readsProcessed++;
				if(maxReads>0 && readsProcessed>maxReads){break;}

				if(Tools.startsWith(line, invalidTitle)){
					valid=false;
					invalidReads++;
					invalidLines++;
					if(deleteInvalid){break;}
				}else{
					assert(Tools.startsWith(line, title)) : "\nline: '"+new String(line)+"'\ntitle: '"+new String(title)+"'\n";
					valid=true;
					validReads++;
					validLines++;
				}
			}else if(gffIn){
				//TODO: Probable bug - after any prior GFF input, later server-mode FASTA sequence lines are classified as valid GFF lines.
				basesProcessed+=line.length;
				valid=!Tools.startsWith(line, invalidGffTitle);
				if(valid){
					validBases+=line.length;
					validLines++;
				}else{
					invalidBases+=line.length;
					invalidLines++;
					//[tax/RenameGiToTaxid#003] FIXED: deleteInvalid gates on invalidReads (processInner_server + process()),
					//but GFF invalid lines previously only bumped invalidLines -> deleteInvalid was a silent no-op for GFF.
					//An unresolvable accession then DELETED the fna (invalidReads>0) but KEPT the gff with tid|-1 seqids ->
					//a broken fna/gff pair that crashes downstream tools assuming matched input. Count it as an invalid read
					//so deleteInvalid drops the gff too, keeping the pair consistent (unresolvable genome skipped entirely).
					invalidReads++;
				}
			}else{
				basesProcessed+=line.length;
				if(valid){
					validBases+=line.length;
					validLines++;
				}else{
					invalidBases+=line.length;
					invalidLines++;
				}
			}

			if(valid || keepAll){
				if(bsw!=null){bsw.println(line);}
			}else{
				if(bswInvalid!=null){bswInvalid.println(line);}
			}
		}
		return valid;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Input file paths to process, preserving command-line order and uniqueness. */
	private LinkedHashSet<String> in1=new LinkedHashSet<String>();
	/** Primary output path for renamed FASTA or GFF records. */
	private String out1=null;
	/** Optional output path for records that could not be assigned a valid taxID. */
	private String outInvalid=null;
	/** Optional output path for unresolved raw headers or GFF seqids. */
	private String badHeaders=null;

	/** Load and query a local taxonomy tree. */
	private boolean useTree=false;
	/** Path to the local GI-to-taxID table. */
	private String giTableFile=null;
	/** Path to local accession-to-taxID input. */
	private String accessionFile=null;
	/** Path to local accession parsing patterns. */
	private String patternFile=null;

	/*--------------------------------------------------------------*/

	/** Maximum number of FASTA records to process, or -1 for unlimited. */
	private long maxReads=-1;

	/** Number of FASTA records with successfully resolved taxonomy information. */
	private long validReads=0;
	/** Total number of bases in valid FASTA records, or bytes in valid GFF lines. */
	private long validBases=0;
	/** Number of FASTA records or GFF features with invalid or missing taxonomy information. */
	private long invalidReads=0;
	/** Total number of bases or bytes in invalid FASTA records or GFF lines. */
	private long invalidBases=0;
	/** Number of unique taxIDs counted when output occurrence suffixes are enabled. */
	private long taxaCounted=0;

	/** Number of nonempty input lines processed in server mode. */
	private long linesIn=0;
	/** Number of GFF or server-mode lines with valid taxonomy information. */
	private long validLines=0;
	/** Number of GFF or server-mode lines with invalid or missing taxonomy information. */
	private long invalidLines=0;

	/** Buffered-byte trigger for server-mode batches; one appended line may exceed it. */
	private long maxStoredBytes=10000000;

	/** FASTA record counter and total FASTA bases or GFF bytes for elapsed-time reporting. */
	private long readsProcessed=0, basesProcessed=0;

	/** Append original header or seqid content after the emitted taxID; required for GFF input. */
	private boolean prefix=true;
	/** Count per-taxID output occurrences instead of preserving the original header when {@link #prefix} is false. */
	private boolean countTable=true;
	/** Retain invalid records in primary output unless an invalid-record output was requested. */
	private boolean keepAll=true;
	/** Truncate FASTA headers at the first SOH character for nr/nt concatenated headers. */
	private boolean shrinkNames=false;
	/** Print unresolved identifier warnings when {@link #badHeaders} is absent. */
	private boolean warnBadHeaders=true;
	/** Query the configured remote tax server instead of local GI or accession tables. */
	private boolean useServer=false;
	/** Maximum unresolved FASTA headers allowed before aborting, or -1 for unlimited. */
	private long maxInvalidHeaders=-1;
	/** Delete the primary output after processing if any invalid record was observed. */
	private boolean deleteInvalid=false;

	/** Current identifier extraction mode for server requests. */
	private int mode;
	/** Server identifier modes: accession token, GI token, whole header, and UNITE token after the first pipe. */
	private static final int ACCESSION_MODE=0, GI_MODE=1, HEADER_MODE=2, UNITE_MODE=3;

	/** Sticky flag set by any GFF input, selecting GFF line/byte reporting and server-batch classification. */
	private boolean gffIn=false;

	/*--------------------------------------------------------------*/

	/** Input file formats matching {@link #in1}. */
	private final ArrayList<FileFormat> ffin1;
	/** Primary output file format. */
	private final FileFormat ffout1;
	/** Optional invalid-record output file format. */
	private final FileFormat ffoutInvalid;
	/** Optional local taxonomy tree used by non-server mode and UNITE server fallbacks. */
	private final TaxTree tree;

	/*--------------------------------------------------------------*/

	/** Destination for progress messages and final statistics. */
	private PrintStream outstream=System.err;
	/** Enables verbose diagnostics in this tool and selected I/O dependencies. */
	public static boolean verbose=false;
	/** True after an input or output stream reports a close error. */
	public boolean errorState=false;
	/** Overwrite existing output files. */
	private boolean overwrite=true;
	/** Append to existing output files instead of overwriting. */
	private boolean append=false;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Mutable header prefix for valid taxonomy assignments. */
	private static byte[] title=">tid|".getBytes();
	/** Header prefix used to mark invalid FASTA records in server mode. */
	private static byte[] invalidTitle=">tid|-1".getBytes();
	/** Seqid prefix used to mark invalid GFF features in server mode. */
	private static byte[] invalidGffTitle="tid|-1".getBytes();
	/** GFF3 sequence-region pragma prefix; its seqid is renamed alongside feature seqids [#004]. */
	private static final String SEQREGION_PRAGMA="##sequence-region";
	/** Exact construction-time rejection for GFF input that would drop feature columns. */
	private static final String GFF_PREFIX_FALSE_MESSAGE=
			"GFF input requires prefix=t; prefix=f would drop feature columns.";

}
