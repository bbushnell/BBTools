package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import align2.QualityTools;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.SamReadStreamer;
import structures.ListNum;
import tracker.ReadStats;
import var2.AnalyzeVars;
import var2.CallVariants;
import var2.ScafMap;
import var2.VarFilter;
import var2.VarMap;
import var2.VcfLoader;

/**
 * @author Brian Bushnell
 * @date Jan 13, 2014
 *
 */
public class CalcTrueQuality {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args){
		ReadStats.COLLECT_QUALITY_STATS=true;
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		CalcTrueQuality x=new CalcTrueQuality(args);
		ReadStats.overwrite=overwrite;
		x.process();
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/** Calls main() but restores original static variable values. */
	public static void main2(String[] args){
		final boolean oldCOLLECT_QUALITY_STATS=ReadStats.COLLECT_QUALITY_STATS;
		final boolean oldoverwrite=ReadStats.overwrite;
		final int oldREAD_BUFFER_LENGTH=Shared.bufferLen();;
		final boolean oldPIGZ=ReadWrite.USE_PIGZ;
		final boolean oldUnPIGZ=ReadWrite.USE_UNPIGZ;
		final int oldZL=ReadWrite.ZIPLEVEL;
		final boolean oldBF1=ByteFile.FORCE_MODE_BF1;
		final boolean oldBF2=ByteFile.FORCE_MODE_BF2;
		final boolean oldTestInterleaved=FASTQ.TEST_INTERLEAVED;
		final boolean oldForceInterleaved=FASTQ.FORCE_INTERLEAVED;
		
		main(args);
		
		ReadStats.COLLECT_QUALITY_STATS=oldCOLLECT_QUALITY_STATS;
		ReadStats.overwrite=oldoverwrite;
		Shared.setBufferLen(oldREAD_BUFFER_LENGTH);
		ReadWrite.USE_PIGZ=oldPIGZ;
		ReadWrite.USE_UNPIGZ=oldUnPIGZ;
		ReadWrite.ZIPLEVEL=oldZL;
		ByteFile.FORCE_MODE_BF1=oldBF1;
		ByteFile.FORCE_MODE_BF2=oldBF2;
		FASTQ.TEST_INTERLEAVED=oldTestInterleaved;
		FASTQ.FORCE_INTERLEAVED=oldForceInterleaved;
	}
	
	public CalcTrueQuality(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
//		SamLine.PARSE_6=false;
//		SamLine.PARSE_7=false;
		SamLine.PARSE_8=false;
//		SamLine.PARSE_10=false;
//		SamLine.PARSE_OPTIONAL=false;
		SamLine.PARSE_OPTIONAL_MD_ONLY=true; //I only need the MD tag..
		
		ReadWrite.USE_PIGZ=false;
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.ZIPLEVEL=8;
//		SamLine.CONVERT_CIGAR_TO_MATCH=true;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("showstats")){
				showStats=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("t") || a.equals("threads")){
				Shared.setThreads(b);
			}else if(a.equals("build") || a.equals("genome")){
				Data.setGenome(Integer.parseInt(b));
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1") || a.equals("sam")){
				assert(b!=null) : "Bad parameter: "+arg;
				in=b.split(",");
			}else if(a.equals("hist") || a.equals("qhist")){
				qhist=b;
			}else if(a.equals("path")){
				Data.setPath(b);
			}else if(a.equals("append") || a.equals("app")){
//				append=ReadStats.append=Parse.parseBoolean(b);
				assert(false) : "This does not work in append mode.";
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("countindels") || a.equals("indels")){
				COUNT_INDELS=Parse.parseBoolean(b);
			}else if(a.equals("callvariants") || a.equals("callvars") || a.equals("callvariations")){
				callVariants=Parse.parseBoolean(b);
			}else if(a.equals("writematrices") || a.equals("write") || a.equals("wm")){
				writeMatrices=Parse.parseBoolean(b);
			}else if(a.equals("ss") || a.equals("samstreamer")){
				if(b!=null && Tools.isDigit(b.charAt(0))){
					useStreamer=true;
					streamerThreads=Tools.max(1, Integer.parseInt(b));
				}else{
					useStreamer=Parse.parseBoolean(b);
				}
			}else if(a.equals("passes") || a.equals("recalpasses")){
				passes=Integer.parseInt(b);
			}
			
			
			else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("ref")){
				ref=b;
				if(ref!=null && !Tools.isReadableFile(ref)) {
					if(ref.equalsIgnoreCase("phix")) {
						ref=Data.findPath("?phix2.fa.gz");
					}
				}
			}else if(a.equals("realign")){
				realign=Parse.parseBoolean(b);
			}else if(a.equals("prefilter")){
				prefilter=Parse.parseBoolean(b);
			}
			
			else if(a.equals("vars") || a.equals("variants") || a.equals("variations") || a.equals("varfile") || a.equals("inv")){
				if(b==null){varFile=null;}
				else if(FileFormat.isVcfFile(a)){vcfFile=b;}
				else{varFile=b;}
			}else if(a.equals("vcf") || a.equals("vcffile")){
				vcfFile=b;
			}else if(filter.parse(a, b, arg)){
				//do nothing
			}
			
			else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQualityAdjust(arg, a, b)){
				//do nothing
			}
			
			else if(in==null && i==0 && Tools.looksLikeInputStream(arg)){
				in=arg.split(",");
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in) || !Tools.testInputFiles(false, true, ref, vcfFile, varFile)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
//			if(ReadWrite.isCompressed(in1)){ByteFile.FORCE_MODE_BF2=true;}
			ByteFile.FORCE_MODE_BF2=true;
		}
		
//		if(!Tools.testOutputFiles(overwrite, append, false, q102out, qbpout, q10out, q12out, qb012out, qb123out, qb234out, qpout, qout, pout)){
//			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+q102out+"\n");
//		}
		threads=Shared.threads();
		if(qhist!=null){readstats=new ReadStats();}
		
		assert(passes==1 || passes==2);
		
		incrQ102=(TRACK_ALL || use_q102[0] || (passes>1 && use_q102[1]));
		incrQap=(TRACK_ALL || use_qap[0] || (passes>1 && use_qap[1]));
		incrQbp=(TRACK_ALL || use_qbp[0] || (passes>1 && use_qbp[1]));
		incrQpt=(TRACK_ALL || use_qpt[0] || (passes>1 && use_qpt[1]));
		incrQbt=(TRACK_ALL || use_qbt[0] || (passes>1 && use_qbt[1]));
		incrQ10=(TRACK_ALL || use_q10[0] || (passes>1 && use_q10[1]));
		incrQ12=(TRACK_ALL || use_q12[0] || (passes>1 && use_q12[1]));
		incrQb12=(TRACK_ALL || use_qb12[0] || (passes>1 && use_qb12[1]));
		incrQb012=(TRACK_ALL || use_qb012[0] || (passes>1 && use_qb012[1]));
		incrQb123=(TRACK_ALL || use_qb123[0] || (passes>1 && use_qb123[1]));
		incrQb234=(TRACK_ALL || use_qb234[0] || (passes>1 && use_qb234[1]));
		incrQ12b12=(TRACK_ALL || use_q12b12[0] || (passes>1 && use_q12b12[1]));
		incrQp=(TRACK_ALL || use_qp[0] || (passes>1 && use_qp[1]));
		incrQ=(TRACK_ALL || use_q[0] || (passes>1 && use_q[1]));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(){
		Timer t=new Timer();
		
		if(callVariants){
			String inString="in="+in[0];
			for(int i=1; i<in.length; i++){
				inString=inString+","+in[i];
			}
			CallVariants cv=new CallVariants(new String[] {inString, "ref="+ref, "realign="+realign, "ploidy="+ploidy, "prefilter="+prefilter});
//			cv.ploidy=ploidy;
//			cv.prefilter=prefilter;
			cv.varFilter.setFrom(filter);
			
			varMap=cv.process(new Timer());
			scafMap=cv.scafMap;
		}else if(varFile!=null || vcfFile!=null){
			if(ref!=null){
				scafMap=ScafMap.loadReference(ref, true);
			}else{
				scafMap=ScafMap.loadSamHeader(in[0]);
			}
			assert(scafMap!=null && scafMap.size()>0) : "No scaffold names were loaded.";
			if(varFile!=null){
				varMap=VcfLoader.loadVarFile(varFile, scafMap);
			}else{
				varMap=VcfLoader.loadVcfFile(vcfFile, scafMap, false, false);
			}
		}
		
		if(varMap==null || varMap.size()==0 || scafMap==null || scafMap.size()==0){
			varMap=null;
			scafMap=null;
		}
		
		SamLine.PARSE_0=true;//Needed for tile
		for(int pass=0; pass<passes; pass++){
			process(pass);
		}
		
		t.stop();
		
		if(showStats){
			readsProcessed/=passes;
			basesProcessed/=passes;
			readsUsed/=passes;
			basesUsed/=passes;
			varsFixed/=passes;
			varsTotal/=passes;

			if(varMap!=null){
				outstream.println(Tools.format("Ignored "+varsFixed+" known variants out of "+varsTotal+" total (%.2f%%).\n", varsFixed*100.0/varsTotal));
			}
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));

			String rpstring=Tools.padKMB(readsUsed, 8);
			String bpstring=Tools.padKMB(basesUsed, 8);

			outstream.println("Reads Used:    "+rpstring);
			outstream.println("Bases Used:    "+bpstring);
		}
		
		if(errorState){
			throw new RuntimeException(this.getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	public void process(final int pass){
		if(pass>0){
			initializeMatrices(pass-1);
		}
		
		int fnum=0;
		for(String s : in){
			process_MT(s, pass, fnum);
			fnum++;
		}
		
		if(writeMatrices){
			writeMatrices(pass);
			gbmatrices.set(pass, null);
		}
		
		System.err.println("Finished pass "+(pass+1)+"\n");
		
		if(errorState){
			throw new RuntimeException(this.getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	public void process_MT(String fname, int pass, int fnum){
		
		assert(gbmatrices.size()==pass || fnum>0) : gbmatrices.size()+", "+pass;
		
		final SamReadStreamer ss;
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff=FileFormat.testInput(fname, FileFormat.SAM, null, true, false);
			if(useStreamer && Shared.threads()>1 && ff.samOrBam()){
				cris=null;
				ss=new SamReadStreamer(ff, streamerThreads, false, maxReads);
				ss.start();
			}else{
				ss=null;
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null);
				if(verbose){System.err.println("Starting cris");}
				cris.start(); //4567
			}
		}
		
		/* Create Workers */
		final int wthreads=Tools.mid(1, threads, 16);
		ArrayList<Worker> alpt=new ArrayList<Worker>(wthreads);
		for(int i=0; i<wthreads; i++){alpt.add(new Worker(cris, ss, pass));}
		for(Worker pt : alpt){pt.start();}
		
		GBMatrixSet gbmatrix;
		if(fnum==0){
			gbmatrix=new GBMatrixSet(pass);
			gbmatrices.add(gbmatrix);
		}else{
			gbmatrix=gbmatrices.get(pass);
			assert(gbmatrix)!=null;
		}
		assert(gbmatrices.size()==pass+1) : gbmatrices.size()+", "+pass;
		
		/* Wait for threads to die, and gather statistics */
		for(int i=0; i<alpt.size(); i++){
			Worker pt=alpt.get(i);
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			alpt.set(i, null);
			
			gbmatrix.add(pt.matrixT);

			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			readsUsed+=pt.readsUsedT;
			basesUsed+=pt.basesUsedT;
			varsFixed+=pt.varsFixedT;
			varsTotal+=pt.varsTotalT;
		}
		
		/* Shut down I/O streams; capture error status */
		errorState|=ReadWrite.closeStreams(cris);
	
	}
	
	static void add(long[] dest, long[] source){
		assert(dest.length==source.length);
		for(int i=0; i<dest.length; i++){dest[i]+=source[i];}
	}
	
	static void add(long[][] dest, long[][] source){
		assert(dest.length==source.length);
		for(int i=0; i<dest.length; i++){add(dest[i], source[i]);}
	}
	
	static void add(long[][][] dest, long[][][] source){
		assert(dest.length==source.length);
		for(int i=0; i<dest.length; i++){add(dest[i], source[i]);}
	}
	
	static void add(long[][][][] dest, long[][][][] source){
		assert(dest.length==source.length);
		for(int i=0; i<dest.length; i++){add(dest[i], source[i]);}
	}
	
	static void add(long[][][][][] dest, long[][][][][] source){
		assert(dest.length==source.length);
		for(int i=0; i<dest.length; i++){add(dest[i], source[i]);}
	}
	
	public void writeMatrices(int pass){
		int oldZL=ReadWrite.ZIPLEVEL;
		ReadWrite.ZIPLEVEL=8;
		gbmatrices.get(pass).write();
		if(qhist!=null){
			readstats=ReadStats.mergeAll();
			readstats.writeQualityToFile(qhist, false);
		}
		ReadWrite.ZIPLEVEL=oldZL;
	}
	
	public static void writeMatrix(String fname, long[][][][][] goodMatrix, long[][][][][] badMatrix, boolean overwrite, boolean append, int pass){
		assert(fname!=null) : "No file specified";
		if(fname.startsWith("?")){
			fname=fname.replaceFirst("\\?", Data.ROOT_QUALITY);
		}
		fname=fname.replace("_p#", "_p"+pass);
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, null, false, overwrite, append, false);
		TextStreamWriter tsw=new TextStreamWriter(ff);
		//System.err.println("Starting tsw for "+fname);
		tsw.start();
		//System.err.println("Started tsw for "+fname);
		StringBuilder sb=new StringBuilder();
		
		final int d0=goodMatrix.length, d1=goodMatrix[0].length, d2=goodMatrix[0][0].length, d3=goodMatrix[0][0][0].length, d4=goodMatrix[0][0][0][0].length;
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
				for(int c=0; c<d2; c++){
					for(int d=0; d<d3; d++){
						for(int e=0; e<d4; e++){
							long good=goodMatrix[a][b][c][d][e];
							long bad=badMatrix[a][b][c][d][e];
							long sum=good+bad;
							if(sum>0){
								sb.append(a);
								sb.append('\t');
								sb.append(b);
								sb.append('\t');
								sb.append(c);
								sb.append('\t');
								sb.append(d);
								sb.append('\t');
								sb.append(e);
								sb.append('\t');
								sb.append(sum);
								sb.append('\t');
								sb.append(bad);
								sb.append('\n');
							}
						}
						if(sb.length()>0){
							tsw.print(sb.toString());
							sb.setLength(0);
						}
					}
				}
			}
		}
		//System.err.println("Writing "+fname);
		tsw.poisonAndWait();
		if(showStats){System.err.println("Wrote "+fname);}
	}
	
	public static void writeMatrix(String fname, long[][][][] goodMatrix, long[][][][] badMatrix, boolean overwrite, boolean append, int pass){
		assert(fname!=null) : "No file specified";
		if(fname.startsWith("?")){
			fname=fname.replaceFirst("\\?", Data.ROOT_QUALITY);
		}
		fname=fname.replace("_p#", "_p"+pass);
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, null, false, overwrite, append, false);
//		assert(false) : new File(fname).canWrite()+", "+new File(fname).getAbsolutePath();
		TextStreamWriter tsw=new TextStreamWriter(ff);
		//System.err.println("Starting tsw for "+fname);
		tsw.start();
		//System.err.println("Started tsw for "+fname);
		StringBuilder sb=new StringBuilder();
		
		final int d0=goodMatrix.length, d1=goodMatrix[0].length, d2=goodMatrix[0][0].length, d3=goodMatrix[0][0][0].length;
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
				for(int c=0; c<d2; c++){
					for(int d=0; d<d3; d++){
						long good=goodMatrix[a][b][c][d];
						long bad=badMatrix[a][b][c][d];
						long sum=good+bad;
						if(sum>0){
							sb.append(a);
							sb.append('\t');
							sb.append(b);
							sb.append('\t');
							sb.append(c);
							sb.append('\t');
							sb.append(d);
							sb.append('\t');
							sb.append(sum);
							sb.append('\t');
							sb.append(bad);
							sb.append('\n');
						}
					}
					if(sb.length()>0){
						tsw.print(sb.toString());
						sb.setLength(0);
					}
				}
			}
		}
		//System.err.println("Writing "+fname);
		tsw.poisonAndWait();
		if(showStats){System.err.println("Wrote "+fname);}
	}
	
	public static void writeMatrix(String fname, long[][][] goodMatrix, long[][][] badMatrix, boolean overwrite, boolean append, int pass){
		assert(fname!=null) : "No file specified";
		if(fname.startsWith("?")){
			fname=fname.replaceFirst("\\?", Data.ROOT_QUALITY);
		}
		fname=fname.replace("_p#", "_p"+pass);
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, null, false, overwrite, append, false);
		TextStreamWriter tsw=new TextStreamWriter(ff);
		//System.err.println("Starting tsw for "+fname);
		tsw.start();
		//System.err.println("Started tsw for "+fname);
		StringBuilder sb=new StringBuilder();
		
		final int d0=goodMatrix.length, d1=goodMatrix[0].length, d2=goodMatrix[0][0].length;
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
				for(int c=0; c<d2; c++){
					long good=goodMatrix[a][b][c];
					long bad=badMatrix[a][b][c];
					long sum=good+bad;
					if(sum>0){
						sb.append(a);
						sb.append('\t');
						sb.append(b);
						sb.append('\t');
						sb.append(c);
						sb.append('\t');
						sb.append(sum);
						sb.append('\t');
						sb.append(bad);
						sb.append('\n');
					}
				}
				if(sb.length()>0){
					tsw.print(sb.toString());
					sb.setLength(0);
				}
			}
		}
		//System.err.println("Writing "+fname);
		tsw.poisonAndWait();
		if(showStats){System.err.println("Wrote "+fname);}
	}
	
	public static void writeMatrix(String fname, long[][] goodMatrix, long[][] badMatrix, boolean overwrite, boolean append, int pass){
		assert(fname!=null) : "No file specified";
		if(fname.startsWith("?")){
			fname=fname.replaceFirst("\\?", Data.ROOT_QUALITY);
		}
		fname=fname.replace("_p#", "_p"+pass);
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, null, false, overwrite, append, false);
		TextStreamWriter tsw=new TextStreamWriter(ff);
		//System.err.println("Starting tsw for "+fname);
		tsw.start();
		//System.err.println("Started tsw for "+fname);
		StringBuilder sb=new StringBuilder();
		
		final int d0=goodMatrix.length, d1=goodMatrix[0].length;
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
				long good=goodMatrix[a][b];
				long bad=badMatrix[a][b];
				long sum=good+bad;
				if(sum>0){
					sb.append(a);
					sb.append('\t');
					sb.append(b);
					sb.append('\t');
					sb.append(sum);
					sb.append('\t');
					sb.append(bad);
					sb.append('\n');
				}
			}
			if(sb.length()>0){
				tsw.print(sb.toString());
				sb.setLength(0);
			}
		}
		//System.err.println("Writing "+fname);
		tsw.poisonAndWait();
		if(showStats){System.err.println("Wrote "+fname);}
	}
	
	public static void writeMatrix(String fname, long[] goodMatrix, long[] badMatrix, boolean overwrite, boolean append, int pass){
		assert(fname!=null) : "No file specified";
		if(fname.startsWith("?")){
			fname=fname.replaceFirst("\\?", Data.ROOT_QUALITY);
		}
		fname=fname.replace("_p#", "_p"+pass);
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, null, false, overwrite, append, false);
		TextStreamWriter tsw=new TextStreamWriter(ff);
		//System.err.println("Starting tsw for "+fname);
		tsw.start();
		//System.err.println("Started tsw for "+fname);
		StringBuilder sb=new StringBuilder();

		final int d0=goodMatrix.length;
		for(int a=0; a<d0; a++){
			long good=goodMatrix[a];
			long bad=badMatrix[a];
			long sum=good+bad;
			if(sum>0){
				sb.append(a);
				sb.append('\t');
				sb.append(sum);
				sb.append('\t');
				sb.append(bad);
				sb.append('\n');
			}
			if(sb.length()>0){
				tsw.print(sb.toString());
				sb.setLength(0);
			}
		}
		//System.err.println("Writing "+fname);
		tsw.poisonAndWait();
		if(showStats){System.err.println("Wrote "+fname);}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int parseTile(String s) {
		//MISEQ08:172:000000000-ABYD0:1:1101:18147:1925
		int colons=0;
		int a=0, b=0;
		for(a=0; colons<4 && a<s.length(); a++) {
			if(s.charAt(a)==':') {colons++;}
		}
		for(b=a; colons<5 && b<s.length(); b++) {
			if(s.charAt(b)==':') {colons++;}
		}
		assert(a<s.length() && s.charAt(a)!=':' && s.charAt(a-1)==':') : a+", "+b+", "+s;
		assert(b<s.length() && s.charAt(b)!=':' && s.charAt(b-1)==':') : a+", "+b+", "+s;
		return Parse.parseInt(s, a, b-1);
	}
	
	public static final void recalibrate(Read r){
		recalibrate(r, true, passes>1);
	}
	
	private static final void recalibrate(Read r, boolean pass0, boolean pass1){
		if(r==null) {return;}
//		System.err.println(r.obj);
//		System.err.println(Arrays.toString(r.quality));
		
		final int pairnum;
		if(USE_PAIRNUM){
			pairnum=r.samline==null ? r.pairnum() : r.samline.pairnum();
		}else{
			pairnum=0;
		}
		final int tile=(USE_TILES ? parseTile(r.name())%TMAX : 0);
		if(pass0){
			byte[] quals2=recalibrate(r.bases, r.quality, pairnum, tile, 0);
			for(int i=0; i<quals2.length; i++){
				r.quality[i]=quals2[i];
			} //Allows calibrating sam output.
		}
		if(pass1){
			byte[] quals2=recalibrate(r.bases, r.quality, pairnum, tile, 1);
			for(int i=0; i<quals2.length; i++){
				r.quality[i]=quals2[i];
			} //Allows calibrating sam output.
		}
		
//		assert(OBSERVATION_CUTOFF==0);
//		assert(false) : pass0+", "+pass1;
//
//		System.err.println(Arrays.toString(r.quality));
//		System.err.println(r.obj);
//		assert(false);
	}
	
	public static final byte[] recalibrate(final byte[] bases, final byte[] quals, final int pairnum, 
			int tile, int pass){
		return cmatrices[pass].recalibrate(bases, quals, pairnum, tile);
	}
	
	public static final void unloadMatrices(){
		for(int i=0; i<passes; i++){
			initialized[i]=false;
			cmatrices[i]=null;
		}
	}
	
	public static final void initializeMatrices(){
		for(int i=0; i<passes; i++){
			initializeMatrices(i);
		}
	}
	
	public static final void initializeMatrices(int pass){
		if(initialized[pass]){return;}
		
		synchronized(initialized){
			if(initialized[pass]){return;}
			assert(cmatrices[pass]==null);
			cmatrices[pass]=new CountMatrixSet(pass);
			cmatrices[pass].load();
			initialized[pass]=true;
		}
		
//		assert(false) : (q102ProbMatrix!=null)+", "+(qbpProbMatrix!=null)+", "+(q10ProbMatrix!=null)+", "+(q12ProbMatrix!=null)+", "+(qb012ProbMatrix!=null)+", "+(qb234ProbMatrix!=null)+", "+(qpProbMatrix!=null);
	}
	
	/*--------------------------------------------------------------*/
	
	private static double modify(final double sum, final double bad, final int phred, final long cutoff){
		double expected=QualityTools.PROB_ERROR[phred];

		double sum2=sum+cutoff;
		double bad2=bad+expected*cutoff;
		double measured=bad2/sum2;

		return measured;
		
//		double modified=Math.pow(measured*measured*measured*expected, 0.25);
////		double modified=Math.sqrt(measured*expected);
////		double modified=(measured+expected)*.5;
//
//		return modified;
	}
	
	public static final float[][][][][] toProbs(long[][][][][] sumMatrix, long[][][][][] badMatrix, final long cutoff){
		final int d0=sumMatrix.length, d1=sumMatrix[0].length, d2=sumMatrix[0][0].length, d3=sumMatrix[0][0][0].length, d4=sumMatrix[0][0][0][0].length;
		float[][][][][] probs=new float[d0][d1][d2][d3][d4];
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
				for(int c=0; c<d2; c++){
					for(int d=0; d<d3; d++){
						for(int e=0; e<d4; e++){
							double sum=sumMatrix[a][b][c][d][e];
							double bad=badMatrix[a][b][c][d][e];
							double modified=modify(sum, bad, b, cutoff);
							probs[a][b][c][d][e]=(float)modified;
						}
					}
				}
			}
		}
		return probs;
	}
	
	public static final float[][][][] toProbs(long[][][][] sumMatrix, long[][][][] badMatrix, final long cutoff){
		final int d0=sumMatrix.length, d1=sumMatrix[0].length, d2=sumMatrix[0][0].length, d3=sumMatrix[0][0][0].length;
		float[][][][] probs=new float[d0][d1][d2][d3];
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
				for(int c=0; c<d2; c++){
					for(int d=0; d<d3; d++){
						double sum=sumMatrix[a][b][c][d];
						double bad=badMatrix[a][b][c][d];
						double modified=modify(sum, bad, b, cutoff);
						probs[a][b][c][d]=(float)modified;
					}
				}
			}
		}
		return probs;
	}
	
	public static final float[][][] toProbs(long[][][] sumMatrix, long[][][] badMatrix, final long cutoff){
		final int d0=sumMatrix.length, d1=sumMatrix[0].length, d2=sumMatrix[0][0].length;
		assert(d0==badMatrix.length) : d0+", "+d1+", "+d2+", "+badMatrix.length;
		assert(d1==badMatrix[0].length) : d0+", "+d1+", "+d2+", "+badMatrix[0].length;
		assert(d2==badMatrix[0][0].length) : d0+", "+d1+", "+d2+", "+badMatrix[0][0].length;
		float[][][] probs=new float[d0][d1][d2];
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
				for(int c=0; c<d2; c++){
					double sum=sumMatrix[a][b][c];
					double bad=badMatrix[a][b][c];
					double modified=modify(sum, bad, b, cutoff);
					probs[a][b][c]=(float)modified;
				}
			}
		}
		return probs;
	}
	
	public static final float[][] toProbs(long[][] sumMatrix, long[][] badMatrix, final long cutoff){
		final int d0=sumMatrix.length, d1=sumMatrix[0].length;
		float[][] probs=new float[d0][d1];
		for(int a=0; a<d0; a++){
			for(int b=0; b<d1; b++){
					double sum=sumMatrix[a][b];
					double bad=badMatrix[a][b];
					double modified=modify(sum, bad, b, cutoff);
					probs[a][b]=(float)modified;
			}
		}
		return probs;
	}
	
	/*--------------------------------------------------------------*/
	
	private static String findPath(String fname){
		assert(fname!=null);
//		return Data.findPath(fname);
		if(fname.startsWith("?")){
			fname=fname.replaceFirst("\\?", Data.ROOT_QUALITY);
		}
		return fname;
	}

	public static final long[][] loadMatrix(String fname, int d0){
		if(fname==null){return null;}
		fname=findPath(fname);
		System.err.println("Loading "+fname+".");

		try{
			long[][] matrix=new long[2][d0];

			TextFile tf=new TextFile(fname, false);
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				String[] split=line.split("\t");
				assert(split.length==3) : Arrays.toString(split);
				int a=Integer.parseInt(split[0]);
				long bases=Long.parseLong(split[1]);
				long errors=Long.parseLong(split[2]);
				matrix[0][a]=bases;
				matrix[1][a]=errors;
			}
			return matrix;
		}catch(RuntimeException e){
			System.err.println("Error - please regenerate calibration matrices.");
			throw(e);
		}
	}

	public static final long[][][] loadMatrix(String fname, int d0, int d1){
		if(fname==null){return null;}
		fname=findPath(fname);
		System.err.println("Loading "+fname+".");

		try{
			long[][][] matrix=new long[2][d0][d1];

			TextFile tf=new TextFile(fname, false);
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				String[] split=line.split("\t");
				assert(split.length==4) : Arrays.toString(split);
				int a=Integer.parseInt(split[0]);
				int b=Integer.parseInt(split[1]);
				long bases=Long.parseLong(split[2]);
				long errors=Long.parseLong(split[3]);
				matrix[0][a][b]=bases;
				matrix[1][a][b]=errors;
			}
			return matrix;
		}catch(RuntimeException e){
			System.err.println("Error - please regenerate calibration matrices.");
			throw(e);
		}
	}

	public static final long[][][][] loadMatrix(String fname, int d0, int d1, int d2){
		if(fname==null){return null;}
		fname=findPath(fname);
		System.err.println("Loading "+fname+".");

		try{
			long[][][][] matrix=new long[2][d0][d1][d2];

			TextFile tf=new TextFile(fname, false);
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				String[] split=line.split("\t");
				assert(split.length==5) : Arrays.toString(split);
				int a=Integer.parseInt(split[0]);
				int b=Integer.parseInt(split[1]);
				int c=Integer.parseInt(split[2]);
				long bases=Long.parseLong(split[3]);
				long errors=Long.parseLong(split[4]);
				matrix[0][a][b][c]=bases;
				matrix[1][a][b][c]=errors;
			}
			return matrix;
		}catch(RuntimeException e){
			System.err.println("Error - please regenerate calibration matrices.");
			throw(e);
		}
	}

	public static final long[][][][][] loadMatrix(String fname, int d0, int d1, int d2, int d3){
		if(fname==null){return null;}
		fname=findPath(fname);
		System.err.println("Loading "+fname+".");

		try{
			long[][][][][] matrix=new long[2][d0][d1][d2][d3];

			TextFile tf=new TextFile(fname, false);
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				String[] split=line.split("\t");
				assert(split.length==6) : Arrays.toString(split);
				int a=Integer.parseInt(split[0]);
				int b=Integer.parseInt(split[1]);
				int c=Integer.parseInt(split[2]);
				int d=Integer.parseInt(split[3]);
				long bases=Long.parseLong(split[4]);
				long errors=Long.parseLong(split[5]);
				matrix[0][a][b][c][d]=bases;
				matrix[1][a][b][c][d]=errors;
			}
			return matrix;
		}catch(RuntimeException e){
			System.err.println("Error - please regenerate calibration matrices.");
			throw(e);
		}
	}

	public static final long[][][][][][] loadMatrix(String fname, int d0, int d1, int d2, int d3, int d4){
		if(fname==null){return null;}
		fname=findPath(fname);
		System.err.println("Loading "+fname+".");

		try{
			long[][][][][][] matrix=new long[2][d0][d1][d2][d3][d4];

			TextFile tf=new TextFile(fname, false);
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				String[] split=line.split("\t");
				assert(split.length==7) : Arrays.toString(split);
				int a=Integer.parseInt(split[0]);
				int b=Integer.parseInt(split[1]);
				int c=Integer.parseInt(split[2]);
				int d=Integer.parseInt(split[3]);
				int e=Integer.parseInt(split[4]);
				long bases=Long.parseLong(split[5]);
				long errors=Long.parseLong(split[6]);
				matrix[0][a][b][c][d][e]=bases;
				matrix[1][a][b][c][d][e]=errors;
			}
			return matrix;
		}catch(RuntimeException e){
			System.err.println("Error - please regenerate calibration matrices.");
			throw(e);
		}
	}
	
	private static byte[] fillBaseToNum(){
		byte[] btn=new byte[128];
		Arrays.fill(btn, (byte)5);
		btn['A']=btn['a']=0;
		btn['C']=btn['c']=1;
		btn['G']=btn['g']=2;
		btn['T']=btn['t']=3;
		btn['U']=btn['u']=3;
		btn['E']=4;
		return btn;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class Worker extends Thread {
		
		Worker(ConcurrentReadInputStream cris_, SamReadStreamer ss_, int pass_){
			cris=cris_;
			ss=ss_;
			pass=pass_;
			matrixT=new GBMatrixSet(pass);
		}
		
		@Override
		public void run(){
			if(cris==null){
				runStreamer();
			}else{
				runCris();
			}
		}
		
		public void runCris(){
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				for(int idx=0; idx<reads.size(); idx++){
					Read r1=reads.get(idx);
					Read r2=r1.mate;
					if(pass>0){
						recalibrate(r1, true, false);
						if(r2!=null){recalibrate(r2, true, false);}
					}
					processLocal(r1);
					processLocal(r2);
				}
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		public void runStreamer(){
			ListNum<Read> ln=ss.nextList();
			ArrayList<Read> reads=(ln==null ? null : ln.list);

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				for(int idx=0; idx<reads.size(); idx++){
					Read r1=reads.get(idx);
					Read r2=r1.mate;
					if(pass>0){
						recalibrate(r1, true, false);
						if(r2!=null){recalibrate(r2, true, false);}
					}
					processLocal(r1);
					processLocal(r2);
				}
				ln=ss.nextList();
				reads=(ln==null ? null : ln.list);
			}
		}
		
//		private void fixVars(Read r, SamLine sl){
//
//			final byte[] match=r.match;
//			final byte[] bases=r.bases;
//			
//			final boolean rcomp=(r.strand()==Shared.MINUS);
//			if(rcomp){r.reverseComplement();}
//			
//			int rpos=sl.pos-1-SamLine.countLeadingClip(sl.cigar, true, true);
//			int scafnum=scafMap.getNumber(sl.rnameS());
//
//			for(int qpos=0, mpos=0; mpos<match.length; mpos++){
//				final byte m=match[mpos];
//				final byte b=bases[qpos];
//				
//				if(m=='S' && scafnum>=0){
//					Var v=new Var(scafnum, rpos, rpos+1, b, Var.SUB);
//					varsTotalT++;
//					if(varMap.containsKey(v)){
//						varsFixedT++;
//						match[mpos]='V';
//					}
//				}
//				
//				if(m!='D'){qpos++;}
//				if(m!='I'){rpos++;}
//			}
//			if(rcomp){r.reverseComplement();}
//		}
		
		private void processLocal(Read r){
			
//			assert(false) : pass+", "+matrixT.pass;
			
			if(r==null){return;}
			final int pairnum;
			final SamLine sl=r.samline;
			if(sl!=null){
				assert(sl.strand()==r.strand());
			}
			if(!USE_PAIRNUM){
				pairnum=0;
			}else if(sl!=null){
				pairnum=sl.pairnum();
			}else{
				pairnum=r.pairnum();
			}
			readsProcessedT++;
			basesProcessedT+=r.length();
			
			if(verbose){outstream.println(r+"\n");}
			
			if(verbose){outstream.println("A");}
			final byte[] quals=r.quality, bases=r.bases;
			if(quals==null || bases==null || r.match==null){return;}
			final boolean needsFixing=(varMap!=null && Read.containsVars(r.match));
			//TODO: Not clear that fixing works right now.

			final int tile=(USE_TILES ? parseTile(r.name())%TMAX : 0);
			
			if(r.shortmatch()){
				r.toLongMatchString(false);
			}
			final byte[] match=r.match;
			if(verbose){outstream.println("B");}
//			if(r.containsNonNMS() || r.containsConsecutiveS(8)){
//				if(verbose){System.err.println("*************************************************** "+new String(match));}
//				return;
//			}
			
			if(needsFixing){
				int x=Read.countVars(r.match, true, true, true);
				int y=AnalyzeVars.fixVars(r, sl, varMap, scafMap);
				varsTotalT+=x;
				varsFixedT+=y;
			}
//			System.err.println(needsFixing);
//			assert(false) : pass+", "+varMap.size();
			
			if(r.strand()==Shared.MINUS){
//				r.reverseComplement();
				Tools.reverseInPlace(match);
			}
			if(verbose){outstream.println("C");}
			
			final byte e='E';
			
			if(readstatsT!=null){
				readstatsT.addToQualityHistogram(r);
			}
			
			readsUsedT++;
			final int aq=Tools.sumInt(quals)/quals.length;
			for(int qpos=0, mpos=0, last=quals.length-1; mpos<match.length; mpos++){
				
				final byte m=match[mpos];
				final byte mprev=match[Tools.max(mpos-1, 0)];
				final byte mnext=match[Tools.min(mpos+1, match.length-1)];
				
//				System.err.println("Processing "+(char)m+"\tqpos="+qpos+"  \tmpos="+mpos);
				
				if(verbose){outstream.print("D");}
				final int q0=(qpos>0 ? Tools.mid(QMAX, quals[qpos-1], 0) : QEND);
				assert(quals!=null && qpos<quals.length) : sl+"\n"+new String(match)+"\n"+(quals==null ? "null" : ""+quals.length)+", "+qpos+", "+match.length+", "+mpos;
				final int q1=quals[qpos];
				final int q2=(qpos<last ? Tools.mid(QMAX, quals[qpos+1], 0) : QEND);
				
				byte b0=qpos>1 ? bases[qpos-2] : e;
				byte b1=qpos>0 ? bases[qpos-1] : e;
				byte b2=bases[qpos];
				byte b3=qpos<last ? bases[qpos+1] : e;
				byte b4=qpos<last-1 ? bases[qpos+2] : e;
				byte n0=baseToNum[b0];
				byte n1=baseToNum[b1];
				byte n2=baseToNum[b2];
				byte n3=baseToNum[b3];
				byte n4=baseToNum[b4];
				
				
				if(m=='N' || !AminoAcid.isFullyDefined(b2)){
					if(verbose){outstream.print("E");}
					//do nothing
				}else if(m=='D'){
					if(verbose){outstream.print("E");}
					//do nothing
				}else if(m=='C'){
					if(verbose){outstream.print("E");}
					//do nothing
				}else{
					final int pos=Tools.min(qpos, LENMAX-1);

					if(verbose){outstream.print("F");}
					basesUsedT++;
					if(m=='m' || (!COUNT_INDELS && m=='I')){
						final int incr;
						if(COUNT_INDELS && (mprev=='D' || mnext=='D')){
							incr=1;
							
							if(incrQ102) matrixT.q102BadMatrix[pairnum][q1][q0][q2]+=1;
							if(incrQap) matrixT.qapBadMatrix[pairnum][q1][aq][pos]+=1;
							if(incrQbp) matrixT.qbpBadMatrix[pairnum][q1][n2][pos]+=1;

							if(incrQpt) matrixT.qptBadMatrix[pairnum][q1][pos][tile]+=1;
							if(incrQbt) matrixT.qbtBadMatrix[pairnum][q1][n2][tile]+=1;

							if(incrQ10) matrixT.q10BadMatrix[pairnum][q1][q0]+=1;
							if(incrQ12) matrixT.q12BadMatrix[pairnum][q1][q0]+=1;
							if(incrQb12) matrixT.qb12BadMatrix[pairnum][q1][n1][n2]+=1;
							if(incrQb012) matrixT.qb012BadMatrix[pairnum][q1][n0][n1][n2]+=1;
							if(incrQb123) matrixT.qb123BadMatrix[pairnum][q1][n1][n2][n3]+=1;
							if(incrQb123) matrixT.qb234BadMatrix[pairnum][q1][n2][n3][n4]+=1;
							if(incrQ12b12) matrixT.q12b12BadMatrix[pairnum][q1][q2][n1][n2]+=1;
							if(incrQp) matrixT.qpBadMatrix[pairnum][q1][pos]+=1;
							if(incrQ) matrixT.qBadMatrix[pairnum][q1]+=1;
							matrixT.pBadMatrix[pairnum][pos]+=1;
							
//							matrixT.q102BadMatrix[pairnum][q1][q0][q2]+=1;
//							matrixT.qbpBadMatrix[pairnum][q1][n2][pos]+=1;
//
//							matrixT.q10BadMatrix[pairnum][q1][q0]+=1;
//							matrixT.q12BadMatrix[pairnum][q1][q0]+=1;
//							matrixT.qb12BadMatrix[pairnum][q1][n1][n2]+=1;
//							matrixT.qb012BadMatrix[pairnum][q1][n0][n1][n2]+=1;
//							matrixT.qb123BadMatrix[pairnum][q1][n1][n2][n3]+=1;
//							matrixT.qb234BadMatrix[pairnum][q1][n2][n3][n4]+=1;
//							matrixT.q12b12BadMatrix[pairnum][q1][q2][n1][n2]+=1;
//							matrixT.qpBadMatrix[pairnum][q1][pos]+=1;
//							matrixT.qBadMatrix[pairnum][q1]+=1;
//							matrixT.pBadMatrix[pairnum][pos]+=1;
						}else{
							incr=2;
						}
						
						if(incrQ102) matrixT.q102GoodMatrix[pairnum][q1][q0][q2]+=incr;
						if(incrQap) matrixT.qapGoodMatrix[pairnum][q1][aq][pos]+=incr;
						
						assert(pairnum<matrixT.qapGoodMatrix.length) : pairnum+", "+matrixT.qapGoodMatrix.length;
						assert(q1<matrixT.qapGoodMatrix[0].length) : q1+", "+matrixT.qapGoodMatrix[0].length;
						assert(n2<matrixT.qapGoodMatrix[0][0].length) : n2+", "+matrixT.qapGoodMatrix[0][0].length;
						assert(pos<matrixT.qapGoodMatrix[0][0][0].length) : pos+", "+matrixT.qapGoodMatrix[0][0][0].length;
						
						if(incrQbp) matrixT.qbpGoodMatrix[pairnum][q1][n2][pos]+=incr;

						if(incrQpt) matrixT.qptGoodMatrix[pairnum][q1][pos][tile]+=1;
						if(incrQbt) matrixT.qbtGoodMatrix[pairnum][q1][n2][tile]+=1;

						if(incrQ10) matrixT.q10GoodMatrix[pairnum][q1][q0]+=incr;
						if(incrQ12) matrixT.q12GoodMatrix[pairnum][q1][q0]+=incr;
						if(incrQb12) matrixT.qb12GoodMatrix[pairnum][q1][n1][n2]+=incr;
						if(incrQb012) matrixT.qb012GoodMatrix[pairnum][q1][n0][n1][n2]+=incr;
						if(incrQb123) matrixT.qb123GoodMatrix[pairnum][q1][n1][n2][n3]+=incr;
						if(incrQb234) matrixT.qb234GoodMatrix[pairnum][q1][n2][n3][n4]+=incr;
						if(incrQ12b12) matrixT.q12b12GoodMatrix[pairnum][q1][q2][n1][n2]+=incr;
						if(incrQp) matrixT.qpGoodMatrix[pairnum][q1][pos]+=incr;
						if(incrQ) matrixT.qGoodMatrix[pairnum][q1]+=incr;
						matrixT.pGoodMatrix[pairnum][pos]+=incr;

//						matrixT.q102GoodMatrix[pairnum][q1][q0][q2]+=incr;
//						matrixT.qbpGoodMatrix[pairnum][q1][n2][pos]+=incr;
//
//						matrixT.q10GoodMatrix[pairnum][q1][q0]+=incr;
//						matrixT.q12GoodMatrix[pairnum][q1][q0]+=incr;
//						matrixT.qb12GoodMatrix[pairnum][q1][n1][n2]+=incr;
//						matrixT.qb012GoodMatrix[pairnum][q1][n0][n1][n2]+=incr;
//						matrixT.qb123GoodMatrix[pairnum][q1][n1][n2][n3]+=incr;
//						matrixT.qb234GoodMatrix[pairnum][q1][n2][n3][n4]+=incr;
//						matrixT.q12b12GoodMatrix[pairnum][q1][q2][n1][n2]+=incr;
//						matrixT.qpGoodMatrix[pairnum][q1][pos]+=incr;
//						matrixT.qGoodMatrix[pairnum][q1]+=incr;
//						matrixT.pGoodMatrix[pairnum][pos]+=incr;
					}else if(m=='S' || m=='I'){
					
//						if(!skip){
							if(incrQ102) matrixT.q102BadMatrix[pairnum][q1][q0][q2]+=2;
							if(incrQap) matrixT.qapBadMatrix[pairnum][q1][aq][pos]+=2;
							if(incrQbp) matrixT.qbpBadMatrix[pairnum][q1][n2][pos]+=2;

							if(incrQpt) matrixT.qptBadMatrix[pairnum][q1][pos][tile]+=2;
							if(incrQbt) matrixT.qbtBadMatrix[pairnum][q1][n2][tile]+=2;

							if(incrQ10) matrixT.q10BadMatrix[pairnum][q1][q0]+=2;
							if(incrQ12) matrixT.q12BadMatrix[pairnum][q1][q0]+=2;
							if(incrQb12) matrixT.qb12BadMatrix[pairnum][q1][n1][n2]+=2;
							if(incrQb012) matrixT.qb012BadMatrix[pairnum][q1][n0][n1][n2]+=2;
							if(incrQb123) matrixT.qb123BadMatrix[pairnum][q1][n1][n2][n3]+=2;
							if(incrQb234) matrixT.qb234BadMatrix[pairnum][q1][n2][n3][n4]+=2;
							if(incrQ12b12) matrixT.q12b12BadMatrix[pairnum][q1][q2][n1][n2]+=2;
							if(incrQp) matrixT.qpBadMatrix[pairnum][q1][pos]+=2;
							if(incrQ) matrixT.qBadMatrix[pairnum][q1]+=2;
							matrixT.pBadMatrix[pairnum][pos]+=2;
//						}

//						matrixT.q102BadMatrix[pairnum][q1][q0][q2]+=2;
//						matrixT.qbpBadMatrix[pairnum][q1][n2][pos]+=2;
//
//						matrixT.q10BadMatrix[pairnum][q1][q0]+=2;
//						matrixT.q12BadMatrix[pairnum][q1][q0]+=2;
//						matrixT.qb12BadMatrix[pairnum][q1][n1][n2]+=2;
//						matrixT.qb012BadMatrix[pairnum][q1][n0][n1][n2]+=2;
//						matrixT.qb123BadMatrix[pairnum][q1][n1][n2][n3]+=2;
//						matrixT.qb234BadMatrix[pairnum][q1][n2][n3][n4]+=2;
//						matrixT.q12b12BadMatrix[pairnum][q1][q2][n1][n2]+=2;
//						matrixT.qpBadMatrix[pairnum][q1][pos]+=2;
//						matrixT.qBadMatrix[pairnum][q1]+=2;
//						matrixT.pBadMatrix[pairnum][pos]+=2;
					}else if(m=='V'){
						match[mpos]='S';
					}else if(m=='i'){
						match[mpos]='I';
					}else if(m=='d'){
						match[mpos]='D';
					}else{
						throw new RuntimeException("Bad symbol m='"+((char)m)+"'\n"+new String(match)+"\n"+new String(bases)+"\n");
					}
				}
				if(m!='D' && m!='d'){qpos++;}
			}
			
		}

		long readsProcessedT=0;
		long basesProcessedT=0;
		final ReadStats readstatsT=(qhist==null ? null : new ReadStats());
		long readsUsedT=0, basesUsedT;
		long varsFixedT=0, varsTotalT=0;
		
		private final ConcurrentReadInputStream cris;
		private final SamReadStreamer ss;
		private final int pass;
		GBMatrixSet matrixT;
		
//		IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
		
	}
	
	/** 
	 * Good Bad Matrix.
	 * Tracks counts of calls being correct or incorrect under the specified conditions.
	 * For example, q12GoodMatrix[0][10][15] tracks the number of times correct calls were observed,
	 * for Q10 bases followed by q15 bases, in read 1 (0 for first dimension).
	 * <br><br>
	 * q = quality<br>
	 * p = position<br>
	 * a = average quality<br>
	 * 1 = this position (other positions are relative to this position)
	 * @author Brian Bushnell
	 * @date Jan 13, 2014
	 *
	 */
	class GBMatrixSet{
		
		GBMatrixSet(int pass_){
			pass=pass_;
			assert(pass==0 || (pass==1));
		}
		
		final void add(GBMatrixSet incr){
			if(incrQ102) CalcTrueQuality.add(q102GoodMatrix, incr.q102GoodMatrix);
			if(incrQap) CalcTrueQuality.add(qapGoodMatrix, incr.qapGoodMatrix);
			if(incrQbp) CalcTrueQuality.add(qbpGoodMatrix, incr.qbpGoodMatrix);
			if(incrQpt) CalcTrueQuality.add(qptGoodMatrix, incr.qptGoodMatrix);
			if(incrQbt) CalcTrueQuality.add(qbtGoodMatrix, incr.qbtGoodMatrix);
			if(incrQ10) CalcTrueQuality.add(q10GoodMatrix, incr.q10GoodMatrix);
			if(incrQ12) CalcTrueQuality.add(q12GoodMatrix, incr.q12GoodMatrix);
			if(incrQb12) CalcTrueQuality.add(qb12GoodMatrix, incr.qb12GoodMatrix);
			if(incrQb012) CalcTrueQuality.add(qb012GoodMatrix, incr.qb012GoodMatrix);
			if(incrQb123) CalcTrueQuality.add(qb123GoodMatrix, incr.qb123GoodMatrix);
			if(incrQb234) CalcTrueQuality.add(qb234GoodMatrix, incr.qb234GoodMatrix);
			if(incrQ12b12) CalcTrueQuality.add(q12b12GoodMatrix, incr.q12b12GoodMatrix);
			if(incrQp) CalcTrueQuality.add(qpGoodMatrix, incr.qpGoodMatrix);
			if(incrQ) CalcTrueQuality.add(qGoodMatrix, incr.qGoodMatrix);
			CalcTrueQuality.add(pGoodMatrix, incr.pGoodMatrix);
			
			if(incrQ102) CalcTrueQuality.add(q102BadMatrix, incr.q102BadMatrix);
			if(incrQap) CalcTrueQuality.add(qapBadMatrix, incr.qapBadMatrix);
			if(incrQbp) CalcTrueQuality.add(qbpBadMatrix, incr.qbpBadMatrix);
			if(incrQpt) CalcTrueQuality.add(qptBadMatrix, incr.qptBadMatrix);
			if(incrQbt) CalcTrueQuality.add(qbtBadMatrix, incr.qbtBadMatrix);
			if(incrQ10) CalcTrueQuality.add(q10BadMatrix, incr.q10BadMatrix);
			if(incrQ12) CalcTrueQuality.add(q12BadMatrix, incr.q12BadMatrix);
			if(incrQb12) CalcTrueQuality.add(qb12BadMatrix, incr.qb12BadMatrix);
			if(incrQb012) CalcTrueQuality.add(qb012BadMatrix, incr.qb012BadMatrix);
			if(incrQb123) CalcTrueQuality.add(qb123BadMatrix, incr.qb123BadMatrix);
			if(incrQb234) CalcTrueQuality.add(qb234BadMatrix, incr.qb234BadMatrix);
			if(incrQ12b12) CalcTrueQuality.add(q12b12BadMatrix, incr.q12b12BadMatrix);
			if(incrQp) CalcTrueQuality.add(qpBadMatrix, incr.qpBadMatrix);
			if(incrQ) CalcTrueQuality.add(qBadMatrix, incr.qBadMatrix);
			CalcTrueQuality.add(pBadMatrix, incr.pBadMatrix);
		}
		
		public void write() {
			if(incrQ102 && q102matrix!=null){writeMatrix(q102matrix, q102GoodMatrix, q102BadMatrix, overwrite, append, pass);}
			if(incrQap && qapmatrix!=null){writeMatrix(qapmatrix, qapGoodMatrix, qapBadMatrix, overwrite, append, pass);}
			if(incrQbp && qbpmatrix!=null){writeMatrix(qbpmatrix, qbpGoodMatrix, qbpBadMatrix, overwrite, append, pass);}
			if(incrQpt && qptmatrix!=null){writeMatrix(qptmatrix, qptGoodMatrix, qptBadMatrix, overwrite, append, pass);}
			if(incrQbt && qbtmatrix!=null){writeMatrix(qbtmatrix, qbtGoodMatrix, qbtBadMatrix, overwrite, append, pass);}
			if(incrQ10 && q10matrix!=null){writeMatrix(q10matrix, q10GoodMatrix, q10BadMatrix, overwrite, append, pass);}
			if(incrQ12 && q12matrix!=null){writeMatrix(q12matrix, q12GoodMatrix, q12BadMatrix, overwrite, append, pass);}
			if(incrQb12 && qb12matrix!=null){writeMatrix(qb12matrix, qb12GoodMatrix, qb12BadMatrix, overwrite, append, pass);}
			if(incrQb012 && qb012matrix!=null){writeMatrix(qb012matrix, qb012GoodMatrix, qb012BadMatrix, overwrite, append, pass);}
			if(incrQb123 && qb123matrix!=null){writeMatrix(qb123matrix, qb123GoodMatrix, qb123BadMatrix, overwrite, append, pass);}
			if(incrQb234 && qb234matrix!=null){writeMatrix(qb234matrix, qb234GoodMatrix, qb234BadMatrix, overwrite, append, pass);}
			if(incrQ12b12 && q12b12matrix!=null){writeMatrix(q12b12matrix, q12b12GoodMatrix, q12b12BadMatrix, overwrite, append, pass);}
			if(incrQp && qpmatrix!=null){writeMatrix(qpmatrix, qpGoodMatrix, qpBadMatrix, overwrite, append, pass);}
			if(incrQ && qmatrix!=null){writeMatrix(qmatrix, qGoodMatrix, qBadMatrix, overwrite, append, pass);}
			if(pmatrix!=null){writeMatrix(pmatrix, pGoodMatrix, pBadMatrix, overwrite, append, pass);}
		}

		final long[][][][] q102GoodMatrix=new long[2][QMAX2][QMAX2][QMAX2];
		final long[][][][] q102BadMatrix=new long[2][QMAX2][QMAX2][QMAX2];

		final long[][][][] qapGoodMatrix=new long[2][QMAX2][QMAX+1][LENMAX];
		final long[][][][] qapBadMatrix=new long[2][QMAX2][QMAX+1][LENMAX];

		final long[][][][] qbpGoodMatrix=new long[2][QMAX2][BMAX][LENMAX];
		final long[][][][] qbpBadMatrix=new long[2][QMAX2][BMAX][LENMAX];

		final long[][][][] qptGoodMatrix=new long[2][QMAX2][LENMAX][TMAX];
		final long[][][][] qptBadMatrix=new long[2][QMAX2][LENMAX][TMAX];
		
		final long[][][][] qbtGoodMatrix=new long[2][QMAX2][BMAX][TMAX];
		final long[][][][] qbtBadMatrix=new long[2][QMAX2][BMAX][TMAX];

		final long[][][] q10GoodMatrix=new long[2][QMAX2][QMAX2];
		final long[][][] q10BadMatrix=new long[2][QMAX2][QMAX2];

		final long[][][] q12GoodMatrix=new long[2][QMAX2][QMAX2];
		final long[][][] q12BadMatrix=new long[2][QMAX2][QMAX2];

		final long[][][][] qb12GoodMatrix=new long[2][QMAX2][BMAX][BMAX];
		final long[][][][] qb12BadMatrix=new long[2][QMAX2][BMAX][BMAX];

		final long[][][][][] qb012GoodMatrix=new long[2][QMAX2][BMAX][BMAX][BMAX];
		final long[][][][][] qb012BadMatrix=new long[2][QMAX2][BMAX][BMAX][BMAX];

		final long[][][][][] qb123GoodMatrix=new long[2][QMAX2][BMAX][BMAX][BMAX];
		final long[][][][][] qb123BadMatrix=new long[2][QMAX2][BMAX][BMAX][BMAX];

		final long[][][][][] qb234GoodMatrix=new long[2][QMAX2][BMAX][BMAX][BMAX];
		final long[][][][][] qb234BadMatrix=new long[2][QMAX2][BMAX][BMAX][BMAX];

		final long[][][][][] q12b12GoodMatrix=new long[2][QMAX2][QMAX2][BMAX][BMAX];
		final long[][][][][] q12b12BadMatrix=new long[2][QMAX2][QMAX2][BMAX][BMAX];

		final long[][][] qpGoodMatrix=new long[2][QMAX2][LENMAX];
		final long[][][] qpBadMatrix=new long[2][QMAX2][LENMAX];

		final long[][] qGoodMatrix=new long[2][QMAX2];
		final long[][] qBadMatrix=new long[2][QMAX2];

		final long[][] pGoodMatrix=new long[2][LENMAX];
		final long[][] pBadMatrix=new long[2][LENMAX];
		
		final int pass;
		
	}
	
	static class CountMatrixSet{
		
		CountMatrixSet(int pass_){
			pass=pass_;
			assert(pass==0 || (pass==1));
			load();
		}
		
		/**
		 * @param bases
		 * @param quals
		 * @param pairnum
		 * @return recalibrated quality scores
		 */
		public byte[] recalibrate(byte[] bases, byte[] quals, int pairnum, int tile) {
			final byte[] quals2;
			final boolean round=(pass<passes-1);
			final int aq=Tools.sumInt(quals)/quals.length;
			if(quals!=null){
				assert(quals.length<=LENMAX || !(use_qp[pass] || use_qbp[pass] || use_qap[pass] || use_qpt[pass])) :
					"\nThese reads are too long ("+quals.length+"bp) for recalibration using position.  Please select different matrices.\n";
				quals2=new byte[quals.length];
				for(int i=0; i<bases.length; i++){
					final byte q2;
					if(!AminoAcid.isFullyDefined(bases[i])){
						q2=0;
					}else{
						final float prob;
						if(USE_SMR){
							prob=estimateErrorProbSMR(quals, bases, i, pairnum, tile, OBSERVATION_CUTOFF[pass], aq);
						}else if(USE_WEIGHTED_AVERAGE){
							prob=estimateErrorProbWeighted(quals, bases, i, pairnum, tile, OBSERVATION_CUTOFF[pass], aq);
						}else if(USE_AVERAGE){
							prob=estimateErrorProbAvg(quals, bases, i, pairnum, aq, tile);
						}else{
							prob=estimateErrorProbMax(quals, bases, i, pairnum, aq, tile);
						}
						q2=Tools.max((byte)2, QualityTools.probErrorToPhred(prob, true));
					}
					quals2[i]=q2;
				}
			}else{
				assert(false) : "Can't recalibrate qualities for reads that don't have quality scores.";
				quals2=null;
				//TODO
			}
			return quals2;
		}

		void load(){
			synchronized(initialized){
				if(initialized[pass]){return;}
				
				if(use_q102[pass]){
					q102CountMatrix=loadMatrix(q102matrix.replace("_p#", "_p"+pass), 2, QMAX2, QMAX2, QMAX2);
					q102ProbMatrix=toProbs(q102CountMatrix[0], q102CountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}
				if(use_qap[pass]){
					qapCountMatrix=loadMatrix(qapmatrix.replace("_p#", "_p"+pass), 2, QMAX2, QMAX+1, LENMAX);
					qapProbMatrix=toProbs(qapCountMatrix[0], qapCountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}
				if(use_qbp[pass]){
					qbpCountMatrix=loadMatrix(qbpmatrix.replace("_p#", "_p"+pass), 2, QMAX2, 4, LENMAX);
					qbpProbMatrix=toProbs(qbpCountMatrix[0], qbpCountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}
				if(use_qpt[pass]){
					qptCountMatrix=loadMatrix(qptmatrix.replace("_p#", "_p"+pass), 2, QMAX2, LENMAX, TMAX);
					qptProbMatrix=toProbs(qptCountMatrix[0], qptCountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}
				if(use_qbt[pass]){
					qbtCountMatrix=loadMatrix(qbtmatrix.replace("_p#", "_p"+pass), 2, QMAX2, 4, TMAX);
					qbtProbMatrix=toProbs(qbtCountMatrix[0], qbtCountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}
				if(use_q10[pass]){
					q10CountMatrix=loadMatrix(q10matrix.replace("_p#", "_p"+pass), 2, QMAX2, QMAX2);
					q10ProbMatrix=toProbs(q10CountMatrix[0], q10CountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}
				if(use_q12[pass]){
					q12CountMatrix=loadMatrix(q12matrix.replace("_p#", "_p"+pass), 2, QMAX2, QMAX2);
					q12ProbMatrix=toProbs(q12CountMatrix[0], q12CountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}
				if(use_qb12[pass]){
					qb12CountMatrix=loadMatrix(qb12matrix.replace("_p#", "_p"+pass), 2, QMAX2, BMAX, 4);
					qb12ProbMatrix=toProbs(qb12CountMatrix[0], qb12CountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}
				if(use_qb012[pass]){
					qb012CountMatrix=loadMatrix(qb012matrix.replace("_p#", "_p"+pass), 2, QMAX2, BMAX, BMAX, 4);
					qb012ProbMatrix=toProbs(qb012CountMatrix[0], qb012CountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}
				if(use_qb123[pass]){
					qb123CountMatrix=loadMatrix(qb123matrix.replace("_p#", "_p"+pass), 2, QMAX2, BMAX, 4, BMAX);
					qb123ProbMatrix=toProbs(qb123CountMatrix[0], qb123CountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}
				if(use_qb234[pass]){
					qb234CountMatrix=loadMatrix(qb234matrix.replace("_p#", "_p"+pass), 2, QMAX2, 4, BMAX, BMAX);
					qb234ProbMatrix=toProbs(qb234CountMatrix[0], qb234CountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}
				if(use_q12b12[pass]){
					q12b12CountMatrix=loadMatrix(q12b12matrix.replace("_p#", "_p"+pass), 2, QMAX2, QMAX2, BMAX, BMAX);
					q12b12ProbMatrix=toProbs(q12b12CountMatrix[0], q12b12CountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}
				if(use_qp[pass]){
					qpCountMatrix=loadMatrix(qpmatrix.replace("_p#", "_p"+pass), 2, QMAX2, LENMAX);
					qpProbMatrix=toProbs(qpCountMatrix[0], qpCountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}
				if(use_q[pass]){
					qCountMatrix=loadMatrix(qmatrix.replace("_p#", "_p"+pass), 2, QMAX2);
					qProbMatrix=toProbs(qCountMatrix[0], qCountMatrix[1], OBSERVATION_CUTOFF[pass]);
				}

				initialized[pass]=true;
			}
		}
		
		public final float estimateErrorProbAvg(byte[] quals, byte[] bases, int pos, int pairnum, int aq, int tile){
			
			final byte e='E';
			final int last=quals.length-1;
			
			final int q0=(pos>0 ? Tools.mid(QMAX, quals[pos-1], 0) : QEND);
			final int q1=quals[pos];
			final int q2=(pos<last ? Tools.mid(QMAX, quals[pos+1], 0) : QEND);
			
			byte b0=pos>1 ? bases[pos-2] : e;
			byte b1=pos>0 ? bases[pos-1] : e;
			byte b2=bases[pos];
			byte b3=pos<last ? bases[pos+1] : e;
			byte b4=pos<last-1 ? bases[pos+2] : e;
			byte n0=baseToNum[b0];
			byte n1=baseToNum[b1];
			byte n2=baseToNum[b2];
			byte n3=baseToNum[b3];
			byte n4=baseToNum[b4];
			
			float expected=PROB_ERROR[q1];
			float sum=0;
			int x=0;

//			System.err.println();
//			System.err.println(((char)b0)+"\t"+((char)b1)+"\t"+((char)b2)+"\t"+((char)b3)+"\t"+((char)b4));
//			System.err.println((n0)+"\t"+(n1)+"\t"+(n2)+"\t"+(n3)+"\t"+(n4));
//			System.err.println(" "+"\t"+(q0)+"\t"+(q1)+"\t"+(q2)+"\t"+(" "));
//			System.err.println("Expected: "+expected);
			
			if(q102ProbMatrix!=null){
				float f=q102ProbMatrix[pairnum][q1][q0][q2];
				sum+=f;
				x++;
			}
			if(qapProbMatrix!=null){
				float f=qapProbMatrix[pairnum][q1][aq][pos];
				sum+=f;
				x++;
			}
			if(qbpProbMatrix!=null){
				float f=qbpProbMatrix[pairnum][q1][n2][pos];
				sum+=f;
				x++;
			}
			if(qptProbMatrix!=null){
				float f=qptProbMatrix[pairnum][q1][pos][tile];
				sum+=f;
				x++;
			}
			if(qbtProbMatrix!=null){
				float f=qbtProbMatrix[pairnum][q1][n2][tile];
				sum+=f;
				x++;
			}
			if(q10ProbMatrix!=null){
				float f=q10ProbMatrix[pairnum][q1][q0];
				sum+=f;
				x++;
			}
			if(q12ProbMatrix!=null){
				float f=q12ProbMatrix[pairnum][q1][q2];
				sum+=f;
				x++;
			}
			if(qb12ProbMatrix!=null){
				float f=qb12ProbMatrix[pairnum][q1][n1][n2];
				sum+=f;
				x++;
			}
			if(qb012ProbMatrix!=null){
				float f=qb012ProbMatrix[pairnum][q1][n0][n1][n2];
				sum+=f;
				x++;
			}
			if(qb123ProbMatrix!=null){
				float f=qb123ProbMatrix[pairnum][q1][n1][n2][n3];
				sum+=f;
				x++;
			}
			if(qb234ProbMatrix!=null){
				float f=qb234ProbMatrix[pairnum][q1][n2][n3][n4];
				sum+=f;
				x++;
			}
			if(q12b12ProbMatrix!=null){
				float f=q12b12ProbMatrix[pairnum][q1][q2][n1][n2];
				sum+=f;
				x++;
			}
			if(qpProbMatrix!=null){
				float f=qpProbMatrix[pairnum][q1][pos];
				sum+=f;
				x++;
			}
			if(qProbMatrix!=null){
				float f=qProbMatrix[pairnum][q1];
				sum+=f;
				x++;
			}
//			System.err.println("result: "+sum+", "+x+", "+sum/(double)x);
//
//			assert(pos<149) : sum+", "+x+", "+sum/(double)x;
			
			if(x<1){
				assert(false);
				return expected;
			}
			return (sum/(float)x);
		}
		
		public final float estimateErrorProbMax(byte[] quals, byte[] bases, int pos, int pairnum, int aq, int tile){
			
			final byte e='E';
			final int last=quals.length-1;
			
			final int q0=(pos>0 ? Tools.mid(QMAX, quals[pos-1], 0) : QEND);
			final int q1=quals[pos];
			final int q2=(pos<last ? Tools.mid(QMAX, quals[pos+1], 0) : QEND);
			
			byte b0=pos>1 ? bases[pos-2] : e;
			byte b1=pos>0 ? bases[pos-1] : e;
			byte b2=bases[pos];
			byte b3=pos<last ? bases[pos+1] : e;
			byte b4=pos<last-1 ? bases[pos+2] : e;
			byte n0=baseToNum[b0];
			byte n1=baseToNum[b1];
			byte n2=baseToNum[b2];
			byte n3=baseToNum[b3];
			byte n4=baseToNum[b4];
			
			final float expected=PROB_ERROR[q1];
			
			float max=-1;
			
			if(q102ProbMatrix!=null){
				float f=q102ProbMatrix[pairnum][q1][q0][q2];
				max=Tools.max(max, f);
			}
			if(qapProbMatrix!=null){
				float f=qapProbMatrix[pairnum][q1][aq][pos];
				max=Tools.max(max, f);
			}
			if(qbpProbMatrix!=null){
				float f=qbpProbMatrix[pairnum][q1][n2][pos];
				max=Tools.max(max, f);
			}
			if(qptProbMatrix!=null){
				float f=qptProbMatrix[pairnum][q1][pos][tile];
				max=Tools.max(max, f);
			}
			if(qbtProbMatrix!=null){
				float f=qbtProbMatrix[pairnum][q1][n2][tile];
				max=Tools.max(max, f);
			}
			if(q10ProbMatrix!=null){
				float f=q10ProbMatrix[pairnum][q1][q0];
				max=Tools.max(max, f);
			}
			if(q12ProbMatrix!=null){
				float f=q12ProbMatrix[pairnum][q1][q2];
				max=Tools.max(max, f);
			}
			if(qb12ProbMatrix!=null){
				float f=qb12ProbMatrix[pairnum][q1][n1][n2];
				max=Tools.max(max, f);
			}
			if(qb012ProbMatrix!=null){
				float f=qb012ProbMatrix[pairnum][q1][n0][n1][n2];
				max=Tools.max(max, f);
			}
			if(qb123ProbMatrix!=null){
				float f=qb123ProbMatrix[pairnum][q1][n1][n2][n3];
				max=Tools.max(max, f);
			}
			if(qb234ProbMatrix!=null){
				float f=qb234ProbMatrix[pairnum][q1][n2][n3][n4];
				max=Tools.max(max, f);
			}
			if(q12b12ProbMatrix!=null){
				float f=q12b12ProbMatrix[pairnum][q1][q2][n1][n2];
				max=Tools.max(max, f);
			}
			if(qpProbMatrix!=null){
				float f=qpProbMatrix[pairnum][q1][pos];
				max=Tools.max(max, f);
			}
			if(qProbMatrix!=null){
				float f=qProbMatrix[pairnum][q1];
				max=Tools.max(max, f);
			}
			
			if(max<0){
				assert(false);
				return expected;
			}
			return max;
		}
		
		public final float estimateErrorProbGeoAvg(byte[] quals, byte[] bases, int pos, int pairnum, int aq, int tile){
			
			final byte e='E';
			final int last=quals.length-1;
			
			final int q0=(pos>0 ? Tools.mid(QMAX, quals[pos-1], 0) : QEND);
			final int q1=quals[pos];
			final int q2=(pos<last ? Tools.mid(QMAX, quals[pos+1], 0) : QEND);
			
			byte b0=pos>1 ? bases[pos-2] : e;
			byte b1=pos>0 ? bases[pos-1] : e;
			byte b2=bases[pos];
			byte b3=pos<last ? bases[pos+1] : e;
			byte b4=pos<last-1 ? bases[pos+2] : e;
			byte n0=baseToNum[b0];
			byte n1=baseToNum[b1];
			byte n2=baseToNum[b2];
			byte n3=baseToNum[b3];
			byte n4=baseToNum[b4];
			
			float expected=PROB_ERROR[q1];
			double product=1;
			int x=0;

//			System.err.println();
//			System.err.println(((char)b0)+"\t"+((char)b1)+"\t"+((char)b2)+"\t"+((char)b3)+"\t"+((char)b4));
//			System.err.println((n0)+"\t"+(n1)+"\t"+(n2)+"\t"+(n3)+"\t"+(n4));
//			System.err.println(" "+"\t"+(q0)+"\t"+(q1)+"\t"+(q2)+"\t"+(" "));
//			System.err.println("Expected: "+expected);
			
			if(q102ProbMatrix!=null){
				float f=q102ProbMatrix[pairnum][q1][q0][q2];
				product*=f;
				x++;
			}
			if(qapProbMatrix!=null){
				float f=qapProbMatrix[pairnum][q1][aq][pos];
				product*=f;
				x++;
			}
			if(qbpProbMatrix!=null){
				float f=qbpProbMatrix[pairnum][q1][n2][pos];
				product*=f;
				x++;
			}
			if(qptProbMatrix!=null){
				float f=qptProbMatrix[pairnum][q1][pos][tile];
				product*=f;
				x++;
			}
			if(qbtProbMatrix!=null){
				float f=qbtProbMatrix[pairnum][q1][n2][tile];
				product*=f;
				x++;
			}
			if(q10ProbMatrix!=null){
				float f=q10ProbMatrix[pairnum][q1][q0];
				product*=f;
				x++;
			}
			if(q12ProbMatrix!=null){
				float f=q12ProbMatrix[pairnum][q1][q2];
				product*=f;
				x++;
			}
			if(qb12ProbMatrix!=null){
				float f=qb12ProbMatrix[pairnum][q1][n1][n2];
				product*=f;
				x++;
			}
			if(qb012ProbMatrix!=null){
				float f=qb012ProbMatrix[pairnum][q1][n0][n1][n2];
				product*=f;
				x++;
			}
			if(qb123ProbMatrix!=null){
				float f=qb123ProbMatrix[pairnum][q1][n1][n2][n3];
				product*=f;
				x++;
			}
			if(qb234ProbMatrix!=null){
				float f=qb234ProbMatrix[pairnum][q1][n2][n3][n4];
				product*=f;
				x++;
			}
			if(q12b12ProbMatrix!=null){
				float f=q12b12ProbMatrix[pairnum][q1][q2][n1][n2];
				product*=f;
				x++;
			}
			if(qpProbMatrix!=null){
				float f=qpProbMatrix[pairnum][q1][pos];
				product*=f;
				x++;
			}
			if(qProbMatrix!=null){
				float f=qProbMatrix[pairnum][q1];
				product*=f;
				x++;
			}
			
			if(x<1){
				assert(false);
				return expected;
			}
			return (float)Math.pow(product, 1.0/x);
		}
		
		public final float estimateErrorProbWeighted(byte[] quals, byte[] bases, int pos, int pairnum, int tile,
				float obs_cutoff, final int aq){
			
			final byte e='E';
			final int last=quals.length-1;
			
			final int q0=(pos>0 ? Tools.mid(QMAX, quals[pos-1], 0) : QEND);
			final int q1=quals[pos];
			final int q2=(pos<last ? Tools.mid(QMAX, quals[pos+1], 0) : QEND);
			
			byte b0=pos>1 ? bases[pos-2] : e;
			byte b1=pos>0 ? bases[pos-1] : e;
			byte b2=bases[pos];
			byte b3=pos<last ? bases[pos+1] : e;
			byte b4=pos<last-1 ? bases[pos+2] : e;
			byte n0=baseToNum[b0];
			byte n1=baseToNum[b1];
			byte n2=baseToNum[b2];
			byte n3=baseToNum[b3];
			byte n4=baseToNum[b4];
			
			long sum=0, bad=0;
			long metrics=0;
			if(q102CountMatrix!=null){
				sum+=q102CountMatrix[0][pairnum][q1][q0][q2];
				bad+=q102CountMatrix[1][pairnum][q1][q0][q2];
				metrics++;
			}
			if(qapCountMatrix!=null){
				sum+=qapCountMatrix[0][pairnum][q1][aq][pos];
				bad+=qapCountMatrix[1][pairnum][q1][aq][pos];
				metrics++;
			}
			if(qbpCountMatrix!=null){
				sum+=qbpCountMatrix[0][pairnum][q1][n2][pos];
				bad+=qbpCountMatrix[1][pairnum][q1][n2][pos];
				metrics++;
			}
			if(qptCountMatrix!=null){
				sum+=qptCountMatrix[0][pairnum][q1][pos][tile];
				bad+=qptCountMatrix[1][pairnum][q1][pos][tile];
				metrics++;
			}
			if(qbtCountMatrix!=null){
				sum+=qbtCountMatrix[0][pairnum][q1][n2][tile];
				bad+=qbtCountMatrix[1][pairnum][q1][n2][tile];
				metrics++;
			}
			if(q10CountMatrix!=null){
				sum+=q10CountMatrix[0][pairnum][q1][q0];
				bad+=q10CountMatrix[1][pairnum][q1][q0];
				metrics++;
			}
			if(q12CountMatrix!=null){
				sum+=q12CountMatrix[0][pairnum][q1][q2];
				bad+=q12CountMatrix[1][pairnum][q1][q2];
				metrics++;
			}
			if(qb12CountMatrix!=null){
				sum+=qb12CountMatrix[0][pairnum][q1][n1][n2];
				bad+=qb12CountMatrix[1][pairnum][q1][n1][n2];
				metrics++;
			}
			if(qb012CountMatrix!=null){
				sum+=qb012CountMatrix[0][pairnum][q1][n0][n1][n2];
				bad+=qb012CountMatrix[1][pairnum][q1][n0][n1][n2];
				metrics++;
			}
			if(qb123CountMatrix!=null){
				sum+=qb123CountMatrix[0][pairnum][q1][n1][n2][n3];
				bad+=qb123CountMatrix[1][pairnum][q1][n1][n2][n3];
				metrics++;
			}
			if(qb234CountMatrix!=null){
				sum+=qb234CountMatrix[0][pairnum][q1][n2][n3][n4];
				bad+=qb234CountMatrix[1][pairnum][q1][n2][n3][n4];
				metrics++;
			}
			if(q12b12CountMatrix!=null){
				sum+=q12b12CountMatrix[0][pairnum][q1][q2][n1][n2];
				bad+=q12b12CountMatrix[1][pairnum][q1][q2][n1][n2];
				metrics++;
			}
			if(qpCountMatrix!=null){
				sum+=qpCountMatrix[0][pairnum][q1][pos];
				bad+=qpCountMatrix[1][pairnum][q1][pos];
				metrics++;
			}
			if(qCountMatrix!=null){
				sum+=qCountMatrix[0][pairnum][q1];
				bad+=qCountMatrix[1][pairnum][q1];
				metrics++;
			}
			
			//TODO: Try taking the sum of roots, then average and square it, instead of the raw numbers.

			final float expectedRate=PROB_ERROR[q1];
			float fakeSum=obs_cutoff;
			float fakeBad=expectedRate*obs_cutoff;
			if(fakeBad<BAD_CUTOFF){
				fakeBad=BAD_CUTOFF;
				fakeSum=BAD_CUTOFF*INV_PROB_ERROR[q1];
			}
			return (float)((bad+fakeBad)/(sum+fakeSum));
		}
		
		public final float estimateErrorProbSMR(byte[] quals, byte[] bases, int pos, 
				int pairnum, int tile, float obs_cutoff, final int aq){
			
			final byte e='E';
			final int last=quals.length-1;
			
			final int q0=(pos>0 ? Tools.mid(QMAX, quals[pos-1], 0) : QEND);
			final int q1=quals[pos];
			final int q2=(pos<last ? Tools.mid(QMAX, quals[pos+1], 0) : QEND);
			
			byte b0=pos>1 ? bases[pos-2] : e;
			byte b1=pos>0 ? bases[pos-1] : e;
			byte b2=bases[pos];
			byte b3=pos<last ? bases[pos+1] : e;
			byte b4=pos<last-1 ? bases[pos+2] : e;
			byte n0=baseToNum[b0];
			byte n1=baseToNum[b1];
			byte n2=baseToNum[b2];
			byte n3=baseToNum[b3];
			byte n4=baseToNum[b4];
			
			double sum=0, bad=0;
			long metrics=0;
			if(q102CountMatrix!=null){
				sum+=Math.sqrt(q102CountMatrix[0][pairnum][q1][q0][q2]);
				bad+=Math.sqrt(q102CountMatrix[1][pairnum][q1][q0][q2]);
				metrics++;
			}
			if(qapCountMatrix!=null){
				sum+=Math.sqrt(qapCountMatrix[0][pairnum][q1][aq][pos]);
				bad+=Math.sqrt(qapCountMatrix[1][pairnum][q1][aq][pos]);
				metrics++;
			}
			if(qbpCountMatrix!=null){
				sum+=Math.sqrt(qbpCountMatrix[0][pairnum][q1][n2][pos]);
				bad+=Math.sqrt(qbpCountMatrix[1][pairnum][q1][n2][pos]);
				metrics++;
			}
			if(qptCountMatrix!=null){
				sum+=Math.sqrt(qptCountMatrix[0][pairnum][q1][pos][tile]);
				bad+=Math.sqrt(qptCountMatrix[1][pairnum][q1][pos][tile]);
				metrics++;
			}
			if(qbtCountMatrix!=null){
				sum+=Math.sqrt(qbtCountMatrix[0][pairnum][q1][n2][tile]);
				bad+=Math.sqrt(qbtCountMatrix[1][pairnum][q1][n2][tile]);
				metrics++;
			}
			if(q10CountMatrix!=null){
				sum+=Math.sqrt(q10CountMatrix[0][pairnum][q1][q0]);
				bad+=Math.sqrt(q10CountMatrix[1][pairnum][q1][q0]);
				metrics++;
			}
			if(q12CountMatrix!=null){
				sum+=Math.sqrt(q12CountMatrix[0][pairnum][q1][q2]);
				bad+=Math.sqrt(q12CountMatrix[1][pairnum][q1][q2]);
				metrics++;
			}
			if(qb12CountMatrix!=null){
				sum+=Math.sqrt(qb12CountMatrix[0][pairnum][q1][n1][n2]);
				bad+=Math.sqrt(qb12CountMatrix[1][pairnum][q1][n1][n2]);
				metrics++;
			}
			if(qb012CountMatrix!=null){
				sum+=Math.sqrt(qb012CountMatrix[0][pairnum][q1][n0][n1][n2]);
				bad+=Math.sqrt(qb012CountMatrix[1][pairnum][q1][n0][n1][n2]);
				metrics++;
			}
			if(qb123CountMatrix!=null){
				sum+=Math.sqrt(qb123CountMatrix[0][pairnum][q1][n1][n2][n3]);
				bad+=Math.sqrt(qb123CountMatrix[1][pairnum][q1][n1][n2][n3]);
				metrics++;
			}
			if(qb234CountMatrix!=null){
				sum+=Math.sqrt(qb234CountMatrix[0][pairnum][q1][n2][n3][n4]);
				bad+=Math.sqrt(qb234CountMatrix[1][pairnum][q1][n2][n3][n4]);
				metrics++;
			}
			if(q12b12CountMatrix!=null){
				sum+=Math.sqrt(q12b12CountMatrix[0][pairnum][q1][q2][n1][n2]);
				bad+=Math.sqrt(q12b12CountMatrix[1][pairnum][q1][q2][n1][n2]);
				metrics++;
			}
			if(qpCountMatrix!=null){
				sum+=Math.sqrt(qpCountMatrix[0][pairnum][q1][pos]);
				bad+=Math.sqrt(qpCountMatrix[1][pairnum][q1][pos]);
				metrics++;
			}
			if(qCountMatrix!=null){
				sum+=Math.sqrt(qCountMatrix[0][pairnum][q1]);
				bad+=Math.sqrt(qCountMatrix[1][pairnum][q1]);
				metrics++;
			}
			
			//TODO: Try taking the sum of roots, then average and square it, instead of the raw numbers.
			
			//Converts back into SquareMeanRoot, probably
			
			//smr=(bad/metrics)^2=bad*bad/(metrics*metrics)=(bad*bad)*(1/(metrics*metrics)
			double mult=metrics<1 ? 1 : 1.0/(metrics*metrics);
			bad=bad*bad*mult;
			sum=sum*sum*mult;
			
			final float expectedRate=PROB_ERROR[q1];
			float fakeSum=obs_cutoff;
			float fakeBad=expectedRate*obs_cutoff;
			if(fakeBad<BAD_CUTOFF){
				fakeBad=BAD_CUTOFF;
				fakeSum=BAD_CUTOFF*INV_PROB_ERROR[q1];
			}
			return (float)((bad+fakeBad)/(sum+fakeSum));
		}
		
		public long[][][][][] q102CountMatrix;
		public long[][][][][] qapCountMatrix;
		public long[][][][][] qbpCountMatrix;
		public long[][][][][] qptCountMatrix;
		public long[][][][][] qbtCountMatrix;
		
		public long[][][][] q10CountMatrix;
		public long[][][][] q12CountMatrix;
		public long[][][][][] qb12CountMatrix;
		public long[][][][][][] qb012CountMatrix;
		public long[][][][][][] qb123CountMatrix;
		public long[][][][][][] qb234CountMatrix;
		public long[][][][][][] q12b12CountMatrix;
		public long[][][][] qpCountMatrix;
		public long[][][] qCountMatrix;

		public float[][][][] q102ProbMatrix;
		public float[][][][] qapProbMatrix;
		public float[][][][] qbpProbMatrix;
		public float[][][][] qptProbMatrix;
		public float[][][][] qbtProbMatrix;
		
		public float[][][] q10ProbMatrix;
		public float[][][] q12ProbMatrix;
		public float[][][][] qb12ProbMatrix;
		public float[][][][][] qb012ProbMatrix;
		public float[][][][][] qb123ProbMatrix;
		public float[][][][][] qb234ProbMatrix;
		public float[][][][][] q12b12ProbMatrix;
		public float[][][] qpProbMatrix;
		public float[][] qProbMatrix;
		
		final int pass;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private VarMap varMap;
	private ScafMap scafMap;
	
	private ReadStats readstats;
	
	private boolean callVariants=false;
	private boolean writeMatrices=true;
	private boolean useStreamer=true;
	private int streamerThreads=3;

	ArrayList<GBMatrixSet> gbmatrices=new ArrayList<GBMatrixSet>();
	
	private PrintStream outstream=System.err;
	private long maxReads=-1;
	private String[] in;
	
	private String qhist=null;
	private String varFile=null;
	private String vcfFile=null;
	
	private long readsProcessed=0;
	private long basesProcessed=0;
	private long readsUsed=0;
	private long basesUsed=0;
	private long varsFixed=0;
	private long varsTotal=0;
	private boolean errorState=false;
	
	private final int threads;
	
	final boolean incrQ102;
	final boolean incrQap;
	final boolean incrQbp;
	final boolean incrQpt;
	final boolean incrQbt;
	final boolean incrQ10;
	final boolean incrQ12;
	final boolean incrQb12;
	final boolean incrQb012;
	final boolean incrQb123;
	final boolean incrQb234;
	final boolean incrQ12b12;
	final boolean incrQp;
	final boolean incrQ;
	
	/*--------------------------------------------------------------*/
	/*----------------         Filter Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	final VarFilter filter=new VarFilter();

	int ploidy=1;
	boolean prefilter=true;
	String ref=null;
	boolean realign=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean showStats=true;
	private static boolean verbose=false;
	private static boolean overwrite=true;
	private static final boolean append=false;
	public static int passes=2;
	
	private static String q102matrix="?q102matrix_p#.txt.gz";
	private static String qapmatrix="?qapmatrix_p#.txt.gz";
	private static String qbpmatrix="?qbpmatrix_p#.txt.gz";
	private static String qptmatrix="?qptmatrix_p#.txt.gz";
	private static String qbtmatrix="?qbtmatrix_p#.txt.gz";
	private static String q10matrix="?q10matrix_p#.txt.gz";
	private static String q12matrix="?q12matrix_p#.txt.gz";
	private static String qb12matrix="?qb12matrix_p#.txt.gz";
	private static String qb012matrix="?qb012matrix_p#.txt.gz";
	private static String qb123matrix="?qb123matrix_p#.txt.gz";
	private static String qb234matrix="?qb234matrix_p#.txt.gz";
	private static String q12b12matrix="?q12b12matrix_p#.txt.gz";
	private static String qpmatrix="?qpmatrix_p#.txt.gz";
	private static String qmatrix="?qmatrix_p#.txt.gz";
	private static String pmatrix="?pmatrix_p#.txt.gz";
	
	private static final boolean[] initialized={false, false};
	
	public static final synchronized void setQmax(int x){
		assert(x>2 && x<94);
		QMAX=x;
		QEND=(QMAX+1);
		QMAX2=(QEND+1);
	}
	private static int QMAX=Read.MAX_CALLED_QUALITY();
	private static int QEND=QMAX+1;
	private static int QMAX2=QEND+1;
	private static final int BMAX=6;
	//Illumina's official specs only go up to 301, I think.
	private static final int LENMAX=361;
	//TMAX=400 works for 10B; 1600 should work for 25B flowcells.
	//10B has 2 swaths per surface and 25B has 4 or 6.  However,
	//1600 makes the matrices too big when there are a lot of threads.
	private static final int TMAX=400;
	private static final byte[] baseToNum=fillBaseToNum();
	private static final byte[] numToBase={'A', 'C', 'G', 'T', 'E', 'N'};
	private static final float[] PROB_ERROR=QualityTools.PROB_ERROR;
	private static final float[] INV_PROB_ERROR=Tools.inverse(PROB_ERROR);
	static{
		PROB_ERROR[0]=0.8f;
		INV_PROB_ERROR[0]=1.25f;
	}
	
	private static final CountMatrixSet[] cmatrices=new CountMatrixSet[2];
	
	public static boolean[] use_q102={false, false};
	public static boolean[] use_qap={false, false};
	public static boolean[] use_qbp={true, true};
	public static boolean[] use_qpt={false, false};
	public static boolean[] use_qbt={false, false};
	public static boolean[] use_q10={false, false};
	public static boolean[] use_q12={false, false};
	public static boolean[] use_qb12={false, false};
	public static boolean[] use_qb012={true, false};
	public static boolean[] use_qb123={true, false};
	public static boolean[] use_qb234={true, false};
	public static boolean[] use_q12b12={false, false};
	public static boolean[] use_qp={false, false};
	public static boolean[] use_q={false, false};
	
	public static boolean USE_SMR=false;
	public static boolean USE_WEIGHTED_AVERAGE=true;
	public static boolean USE_AVERAGE=true;
	public static boolean USE_PAIRNUM=true;
	public static boolean COUNT_INDELS=true;
	public static boolean TRACK_ALL=false;
	public static boolean USE_TILES=use_qpt[0] || use_qpt[1] || use_qbt[0] || use_qbt[1];
	
	public static long OBSERVATION_CUTOFF[]={100, 200}; //Soft threshold
	public static float BAD_CUTOFF=0.5f; //Soft threshold
	
}
