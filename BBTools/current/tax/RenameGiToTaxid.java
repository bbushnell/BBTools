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
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
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
 * @author Brian Bushnell
 * @date Mar 10, 2015
 *
 */
public class RenameGiToTaxid {
	
	public static void main(String[] args){
		Timer t=new Timer();
		RenameGiToTaxid x=new RenameGiToTaxid(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
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
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
			}else if(a.equals("accession")){
				accessionFile=b;
			}else if(a.equals("pattern")){
				patternFile=b;
			}else if(a.equals("tree") || a.equals("taxtree")){
				taxTreeFile=b;
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
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(useServer){
			giTableFile=null;
			accessionFile=null;
			patternFile=null;
			if(mode!=UNITE_MODE){taxTreeFile=null;}
		}//else if taxpath!=null... set them
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out1=parser.out1;
		}
		
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
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
		
		if(ffoutInvalid!=null){keepAll=false;}
		
		assert(giTableFile!=null || accessionFile!=null || TaxTree.SILVA_MODE || useServer) : "No gi or accession information loaded.";
		
		if(taxTreeFile!=null){
			tree=TaxTree.loadTaxTree(taxTreeFile, outstream, true, false);
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
//			System.gc();
		}
	}
	
	void process(Timer t){
		
		ByteStreamWriter bsw=(ffout1==null ? null : new ByteStreamWriter(ffout1)); //Actually, this is required.
		if(bsw!=null){bsw.start();}
		
		ByteStreamWriter bswInvalid=null;
		if(ffoutInvalid!=null){
			bswInvalid=new ByteStreamWriter(ffoutInvalid);
			bswInvalid.start();
		}
		
		ByteStreamWriter bswBadHeaders=null;
		if(badHeaders!=null) {
			bswBadHeaders=new ByteStreamWriter(badHeaders, overwrite, append, false);
			bswBadHeaders.start();
		}
		
		final HashArray1D counts=(countTable && !prefix) ? new HashArray1D(256000, -1L, true) : null;
		
		gffIn=false;
		for(FileFormat ffin : ffin1){
			gffIn=gffIn||ffin.gff();
			ByteFile bf=ByteFile.makeByteFile(ffin);
			if(useServer){
				processInner_server(bf, bsw, bswInvalid, bswBadHeaders, counts, ffin.format());
			}else{
//				IntList list=(useServer ? getIds(bf) : null);
				processInner(bf, bsw, bswInvalid, bswBadHeaders, counts, null);
			}
		}
		
		if(bsw!=null){
			errorState|=bsw.poisonAndWait();
			if(deleteInvalid && invalidReads>0 && !ffout1.stdio()){
				try {
					System.err.println("Deleting "+out1);
					new File(out1).delete();
				} catch (Exception e) {
					System.err.println("An error occured while attempting to delete "+out1);
					e.printStackTrace();
				}
			}
		}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}
		if(bswBadHeaders!=null){errorState|=bswBadHeaders.poisonAndWait();}
		
		t.stop();
		if(!gffIn) {
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
	
	//Unused; not efficient
//	public IntList getIds(ByteFile bf){
//		IntList ids=new IntList();
//		
//		int readsProcessedInner=0;
//		
//		byte[] line=bf.nextLine();
//		ByteBuilder bb=new ByteBuilder();
//		while(line!=null){
//			if(line.length>0 && line[0]=='>'){
//				readsProcessedInner++;
//				if(maxReads>0 && readsProcessedInner>maxReads){break;}
//				
//				for(int i=1; i<line.length; i++){
//					byte b=line[i];
//					if(b==' ' || b=='.'){break;}
//					else{bb.append(b);}
//				}
//				bb.append(',');
//				if(bb.length()>100000){
//					bb.setLength(bb.length()-1);
//					int[] ret;
//					if(mode==ACCESSION_MODE){
//						ret=TaxClient.accessionToTaxidArray(bb.toString());
//					}else if(mode==GI_MODE){
//						ret=TaxClient.giToTaxidArray(bb.toString());
//					}else{
//						ret=TaxClient.headerToTaxidArray(bb.toString());
//					}
//					assert(ret!=null) : bb.toString();
//					for(int i : ret){ids.add(i);}
//					bb.clear();
//				}
//			}
//			line=bf.nextLine();
//		}
//		if(bb.length()>0){
//			bb.setLength(bb.length()-1);
//			int[] ret;
//			if(mode==ACCESSION_MODE){
//				ret=TaxClient.accessionToTaxidArray(bb.toString());
//			}else if(mode==GI_MODE){
//				ret=TaxClient.giToTaxidArray(bb.toString());
//			}else{
//				ret=TaxClient.headerToTaxidArray(bb.toString());
//			}
//			assert(ret!=null) : bb.toString();
//			for(int i : ret){ids.add(i);}
//			bb.clear();
//		}
//		
//		bf.reset();
//		return ids;
//	}
	
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
						
						if(shrinkNames){//This is for nr/nt
							for(int i=initial; i<terminal; i++){
								if(line[i]==1){//SOH
									terminal=i;
								}
							}
						}

						String s=new String(line, initial, terminal-initial);

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
	
	private void updateHeadersFromServer(ArrayList<byte[]> lines, HashArray1D counts, ByteStreamWriter bswBadHeaders, int format){
		if(format==FileFormat.FA){
			updateHeadersFromServer_fasta(lines, counts, bswBadHeaders);
		}else if(format==FileFormat.GFF){
			updateHeadersFromServer_gff(lines, counts, bswBadHeaders);
		}else{
			assert(false) : "Unsupported type: "+format;
		}
	}
	
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
		
//		System.err.println("Sending '"+bb+"'");
		
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
//					else {
//						assert(false) : names.get(i);
//					}
				}else{
					//Sometimes the species gets renamed.
//					TaxNode tn=tree.getNodeByName(name);
//					if(tn==null || tn.id==serverIds[i]) {System.err.println(name+", "+serverIds[i]+", "+tn+", "+tree.getNodesByName(name));}
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
			KillSwitch.kill("Maximum bad headers exceeded: "+maxInvalidHeaders);
		}
	}
	
	private int[] translateToTaxIDs(ByteBuilder bb, int terms) {
		//This block was created due to mysterious crashes processing nr.faa.
		//They seemed to be due to incomplete responses.
		int[] response=null;
		final String query=bb.toString();
		for(int i=0; i<10 && response==null; i++) {
			if(mode==ACCESSION_MODE || mode==UNITE_MODE){
				response=TaxClient.accessionToTaxidArray(query);
			}else if(mode==GI_MODE){
				response=TaxClient.giToTaxidArray(query);
			}else{
				response=TaxClient.headerToTaxidArray(query);//TODO: Keep this as a ByteBuilder
			}
			if(response==null) {
				System.err.println("Retrying null response for length "+terms+" query: "+query);
			}else if(response.length!=terms) {
				System.err.println("Retrying corrupt length "+response.length+" response for length "+terms+" query: "+query);
				response=null;
			}
			if(response==null) {
				try {
					Thread.sleep(500*(i+2));
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
		return response;
	}
	
	private void updateHeadersFromServer_gff(ArrayList<byte[]> lines, HashArray1D counts, ByteStreamWriter bswBadHeaders){
		ByteBuilder bb=new ByteBuilder();
		ArrayList<String> names=new ArrayList<String>();
		for(byte[] line : lines){
			if(line[0]!='#' && !Tools.startsWith(line, "tid")){
				if(bb.length()>0){bb.append(',');}
				for(byte b : line){
					if(b=='\t'){break;}
					bb.append(b);
				}
			}
		}
		if(bb.length()<1){return;}
		
//		assert(false) : bb;
		
//		System.err.println("Sending '"+bb+"'");
		
		int[] serverIds;
		if(mode==ACCESSION_MODE || mode==UNITE_MODE){
			serverIds=TaxClient.accessionToTaxidArray(bb.toString());
		}else if(mode==GI_MODE){
			serverIds=TaxClient.giToTaxidArray(bb.toString());
		}else{
			serverIds=TaxClient.headerToTaxidArray(bb.toString());
		}
		if(serverIds==null){
			KillSwitch.kill("Null response for '"+bb.toString()+"'");
		}
//		assert(serverIds!=null) : "Null response for '"+bb.toString()+"'";
		bb.clear();
		
		if(!names.isEmpty()){
			assert(tree!=null) : "Need to load a TaxTree.";
			assert(names.size()==serverIds.length);
			for(int i=0; i<serverIds.length; i++){
				final String name=names.get(i);
				if(serverIds[i]<0){
					TaxNode tn=tree.getNodeByName(name);
					if(tn!=null){serverIds[i]=tn.id;}
//					else {
//						assert(false) : names.get(i);
//					}
				}else{
					//Sometimes the species gets renamed.
//					TaxNode tn=tree.getNodeByName(name);
//					if(tn==null || tn.id==serverIds[i]) {System.err.println(name+", "+serverIds[i]+", "+tn+", "+tree.getNodesByName(name));}
				}
			}
		}
		
		for(int lineNum=0, serverNum=0; lineNum<=lines.size(); lineNum++){
			byte[] line=lines.get(lineNum);
			if(line[0]!='#' && !Tools.startsWith(line, "tid")){
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
				if(serverNum>=serverIds.length){break;}
			}
		}
		if(maxInvalidHeaders>=0 && invalidReads>maxInvalidHeaders){
			KillSwitch.kill("Maximum bad headers exceeded: "+maxInvalidHeaders);
		}
	}
	
	private void processInner_server(ByteFile bf, ByteStreamWriter bsw, ByteStreamWriter bswInvalid, ByteStreamWriter bswBadHeaders, HashArray1D counts, int format){
		
		ArrayList<byte[]> lines=new ArrayList<byte[]>();
		byte[] line=bf.nextLine();
		boolean valid=false;
		long storedBytes=0;
		
		while(line!=null){
			
			if(line.length>0){
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
					//Somehow this line is firing:
//					Exception in thread "main" java.lang.AssertionError:
//						line: '>WP_330397346.1 spore protease YyaC [Maledivibacter halophilus]SKC56058.1 putative sporulation protein YyaC [Maledivibacter halophilus]'
//						title: '>tid|'
//
//						        at tax.RenameGiToTaxid.dumpBuffer(RenameGiToTaxid.java:790)
//						        at tax.RenameGiToTaxid.processInner_server(RenameGiToTaxid.java:751)
//						        at tax.RenameGiToTaxid.process(RenameGiToTaxid.java:262)
//						        at tax.RenameGiToTaxid.main(RenameGiToTaxid.java:39)
					assert(Tools.startsWith(line, title)) : "\nline: '"+new String(line)+"'\ntitle: '"+new String(title)+"'\n";
					valid=true;
					validReads++;
					validLines++;
				}
			}else if(gffIn){
				basesProcessed+=line.length;
				valid=!Tools.startsWith(line, invalidGffTitle);
				if(valid){
					validBases+=line.length;
					validLines++;
				}else{
					invalidBases+=line.length;
					invalidLines++;
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
	
	
	/*--------------------------------------------------------------*/
	
	private LinkedHashSet<String> in1=new LinkedHashSet<String>();
	private String out1=null;
	private String outInvalid=null;
	private String badHeaders=null;

	private String taxTreeFile=null;
	private String giTableFile=null;
	private String accessionFile=null;
	private String patternFile=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;

	private long validReads=0;
	private long validBases=0;
	private long invalidReads=0;
	private long invalidBases=0;
	private long taxaCounted=0;

	private long linesIn=0;
	private long validLines=0;
	private long invalidLines=0;
	
	private long maxStoredBytes=10000000;
	
	private long readsProcessed=0, basesProcessed=0;

	private boolean prefix=true;
	private boolean countTable=true;
	private boolean keepAll=true;
	private boolean shrinkNames=false;
	private boolean warnBadHeaders=true;
	private boolean useServer=false;
	/** Crash if the number of invalid headers exceeds this */
	private long maxInvalidHeaders=-1;
	/** Delete the output file if there are any invalid headers */
	private boolean deleteInvalid=false;
	
	private int mode;
	private static final int ACCESSION_MODE=0, GI_MODE=1, HEADER_MODE=2, UNITE_MODE=3;
	
	private boolean gffIn=false;
	
	/*--------------------------------------------------------------*/
	
	private final ArrayList<FileFormat> ffin1;
	private final FileFormat ffout1;
	private final FileFormat ffoutInvalid;
	private final TaxTree tree;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;

	private static byte[] title=">tid|".getBytes();
	private static byte[] invalidTitle=">tid|-1".getBytes();
	private static byte[] invalidGffTitle="tid|-1".getBytes();
	
}
