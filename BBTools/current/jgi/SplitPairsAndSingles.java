package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.DualCris;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import structures.ListNum;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Sep 4, 2013
 *
 */
public final class SplitPairsAndSingles {
	
	public static void main(String[] args){
		SplitPairsAndSingles x=new SplitPairsAndSingles(args);
		x.process();
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	public SplitPairsAndSingles(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		Parser parser=new Parser();
		parser.trimq=trimq;
		parser.minReadLength=minReadLength;
		boolean setOut=false, setOuts=false;
		boolean fixInterleaving_=false, repair_=false, allowIdenticalPairNames_=false;

		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("out") || a.equals("out1") || a.equals("outp") || a.equals("outp1") || a.equals("outpair") || a.equals("outpair1")){
				out1=b;
				setOut=true;
			}else if(a.equals("out2") || a.equals("outp2") || a.equals("outpair2")){
				out2=b;
			}else if(a.equals("outs") || a.equals("outsingle") || a.equals("outb") || a.equals("outbad")){
				outsingle=b;
				setOut=true;
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("addslash")){
				addslash=Parse.parseBoolean(b);
			}else if(a.equals("addcolon")){
				addcolon=Parse.parseBoolean(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("fixinterleaving") || a.equals("fi") || a.equals("fint") || a.equals("fixint")){
				fixInterleaving_=Parse.parseBoolean(b);
				if(fixInterleaving_){repair_=false;}
			}else if(a.equals("allowidenticalnames") || a.equals("ain")){
				allowIdenticalPairNames_=Parse.parseBoolean(b);
			}else if(a.equals("repair") || a.equals("rp")){
				repair_=Parse.parseBoolean(b);
				if(repair_){fixInterleaving_=false;}
			}else if(i==0 && in1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				in1=args[i];
			}else if(i==1 && out1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				out1=args[i];
				setOut=true;
			}else if(i==2 && outsingle==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				outsingle=args[i];
				setOuts=true;
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=parser.trimq;
			trimE=parser.trimE();
			minReadLength=parser.minReadLength;
		}

		allowIdenticalPairNames=allowIdenticalPairNames_;
		fixInterleaving=fixInterleaving_;
		repair=repair_;
		assert(!repair || ! fixInterleaving) : "ERROR: Choose 'fixInterleaving' or 'repair', but not both.";
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(in1!=null && in1.contains("#") && !new File(in1).exists()){
			int pound=in1.lastIndexOf('#');
			String a=in1.substring(0, pound);
			String b=in1.substring(pound+1);
			in1=a+1+b;
			in2=a+2+b;
		}
		if(in2!=null && (in2.contains("=") || in2.equalsIgnoreCase("null"))){in2=null;}
		
		if(fixInterleaving){
			if(in2!=null){
				System.err.println("ERROR: 'FixInterleaving' mode only works with a single interleaved input file, not paired input files.");
				System.err.println("Aborting.");
				System.exit(1);
			}
			parser.setInterleaved=true;
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
			outstream.println("Paired input disabled; running in FixInterleaving mode");
		}
		
		if(repair){
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
			outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
		}else{
			if(!parser.setInterleaved && in2==null){
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=true;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}
			if(in2!=null){
				if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
			}
		}
		
		if(out1!=null && out1.contains("#")){
			int pound=out1.lastIndexOf('#');
			String a=out1.substring(0, pound);
			String b=out1.substring(pound+1);
			out1=a+1+b;
			out2=a+2+b;
		}

		if(!setOut){
			System.err.println("No output stream specified.  To write to stdout, please specify 'out=stdout.fq' or similar.");
//			out1="stdout.fq";
			outstream=System.err;
			out2=null;
		}else if("stdout".equalsIgnoreCase(out1) || "standarddout".equalsIgnoreCase(out1)){
			out1="stdout.fq";
			outstream=System.err;
			out2=null;
		}
		if(out1!=null && !Tools.canWrite(out1, overwrite)){throw new RuntimeException("Output file "+out1+" already exists, and overwrite="+overwrite);}

		assert(!in1.equalsIgnoreCase(out1));
		assert(!in1.equalsIgnoreCase(outsingle));
		assert(!in1.equalsIgnoreCase(in2));
		assert(out1==null || !out1.equalsIgnoreCase(out2)) : "out2 may not be defined without out1, and out1 may not equal out2.";
		assert(out1==null || !out1.equalsIgnoreCase(outsingle));
		
		pairMap=(repair ? new LinkedHashMap<String, Read>() : null);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}

	public void process(){
		
		Timer t=new Timer();
		
		process2();
		
		t.stop();
		
		outstream.println("\nInput:                  \t"+readsIn+" reads \t\t"+basesIn+" bases.");
		
		if(qtrimLeft || qtrimRight){
			outstream.println("Trimmed:                \t"+readsTrimmed+" reads ("+Tools.format("%.2f",readsTrimmed*100.0/readsIn)+"%) \t"+
					basesTrimmed+" bases ("+Tools.format("%.2f",basesTrimmed*100.0/basesIn)+"%)");
		}
		outstream.println("Result:                 \t"+readsOut+" reads ("+Tools.format("%.2f",readsOut*100.0/readsIn)+"%) \t"+
				basesOut+" bases ("+Tools.format("%.2f",basesOut*100.0/basesIn)+"%)");
		outstream.println("Pairs:                  \t"+pairsOut+" reads ("+Tools.format("%.2f",pairsOut*100.0/readsIn)+"%) \t"+
				pairBasesOut+" bases ("+Tools.format("%.2f",pairBasesOut*100.0/basesIn)+"%)");
		outstream.println("Singletons:             \t"+singlesOut+" reads ("+Tools.format("%.2f",singlesOut*100.0/readsIn)+"%) \t"+
				singleBasesOut+" bases ("+Tools.format("%.2f",singleBasesOut*100.0/basesIn)+"%)");
		
		if(showSpeed){
			outstream.println();
			outstream.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private void process2(){
		final ConcurrentReadInputStream cris;
		if(in2!=null && repair){
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=DualCris.getReadInputStream(maxReads, true, ff1, ff2, null, null);
		}else{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, null, null, null);
		}
		if(verbose){System.err.println("Started cris");}
		cris.start();
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
		final ConcurrentReadOutputStream ros, rosb;
		final int buff=4;
		if(out1!=null){
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, false);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, false);
			ros=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
			ros.start();
		}else{ros=null;}
		if(outsingle!=null){
			FileFormat ff1=FileFormat.testOutput(outsingle, FileFormat.FASTQ, null, true, overwrite, append, false);
			rosb=ConcurrentReadOutputStream.getStream(ff1, null, buff, null, true);
			rosb.start();
		}else{rosb=null;}
		if(ros!=null || rosb!=null){
			outstream.println("Started output stream.");
		}
		
//		assert(false) : out1+", "+out2+", "+outsingle;
		if(fixInterleaving){
			process3_fixInterleaving(cris, ros, rosb);
		}else if(repair){
			if(cris.getClass()==DualCris.class){
				process3_repair((DualCris)cris, ros, rosb);
			}else{
				process3_repair(cris, ros, rosb);
			}
		}else{
			process3(cris, ros, rosb);
		}

		
		ReadWrite.closeStreams(cris, ros, rosb);
	}
//
//	private void process3_old(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros, final ConcurrentReadOutputStream rosb){
//
//		ListNum<Read> ln=cris.nextList();
//		ArrayList<Read> reads0=(ln!=null ? ln.list : null);
//		ArrayList<Read> single=(rosb==null ? null : new ArrayList<Read>(Shared.bufferLen()));
//
//		while(reads0!=null && reads0.size()>0){
//			ArrayList<Read> reads=(ArrayList<Read>) reads0.clone();
//			int removed=0;
//			for(int i=0; i<reads.size(); i++){
//				Read r1=reads.get(i);
//				Read r2=r1.mate;
//
//				readsIn++;
//				basesIn+=r1.length();
//				if(r2!=null){
//					readsIn++;
//					basesIn+=r2.length();
//				}
//
//				{
//					if(trimLeft || trimRight){
//						if(r1!=null){
//							int x=TrimRead.trimFast(r1, trimLeft, trimRight, trimq, 1);
//							basesTrimmed+=x;
//							readsTrimmed+=(x>0 ? 1 : 0);
//						}
//						if(r2!=null){
//							int x=TrimRead.trimFast(r2, trimLeft, trimRight, trimq, 1);
//							basesTrimmed+=x;
//							readsTrimmed+=(x>0 ? 1 : 0);
//						}
//					}
//
//					final int rlen1=(r1==null ? -1 : r1.length());
//					final int rlen2=(r2==null ? -1 : r2.length());
//
//					if(verbose){System.err.println("rlen1="+rlen1+", rlen2="+rlen2);}
//
//					if(rlen1<minReadLength || rlen2<minReadLength){
//						reads.set(i, null);
//						removed++;
//						r1.mate=null;
//						if(r2!=null){
//							r2.mate=null;
//						}
//
//						if(rlen1>=minReadLength){
//							single.add(r1);
//							singlesOut++;
//							singleBasesOut+=rlen1;
//						}
//						if(rlen2>=minReadLength){
//							single.add(r2);
//							singlesOut++;
//							singleBasesOut+=rlen2;
//						}
//					}else{
//						if(r1!=null){
//							pairsOut++;
//							pairBasesOut+=rlen2;
//						}
//						if(r2!=null){
//							pairsOut++;
//							pairBasesOut+=rlen2;
//						}
//					}
//				}
//			}
//
//			if(rosb!=null){
//				if(verbose){System.err.println("Adding "+single.size()+" to single out.");}
//				rosb.add(new ArrayList<Read>(single), ln.id);
//				single.clear();
//			}
//
//			if(ros!=null){
//				if(removed>0){Tools.condenseStrict(reads);}
//				ArrayList<Read> x=new ArrayList<Read>(reads.size());
//				x.addAll(reads);
//				if(verbose){System.err.println("Adding "+x.size()+" to pair out.");}
//				ros.add(x, ln.id);
//			}
//
//			cris.returnList(ln);
//			ln=cris.nextList();
//			reads0=(ln!=null ? ln.list : null);
//		}
//		cris.returnList(ln);
//
//		readsOut+=singlesOut+pairsOut;
//		basesOut+=singleBasesOut+pairBasesOut;
//	}
	
	private void process3(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros, final ConcurrentReadOutputStream rosb){

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=ln.list;
		
		final ArrayList<Read> pairs=(ros==null ? null : new ArrayList<Read>(Shared.bufferLen()));
		final ArrayList<Read> singles=(rosb==null ? null : new ArrayList<Read>(Shared.bufferLen()));
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			for(int i=0; i<reads.size(); i++){
				Read r1=reads.get(i);
				Read r2=r1.mate;
				processPair(r1, r2, pairs, singles);
			}
			
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
			
			if(rosb!=null){
				if(verbose){System.err.println("Adding "+singles.size()+" to single out.");}
				rosb.add(new ArrayList<Read>(singles), ln.id);
				singles.clear();
			}
			
			if(ros!=null){
				if(verbose){System.err.println("Adding "+pairs.size()+" to pair out.");}
				ros.add(new ArrayList<Read>(pairs), ln.id);
				pairs.clear();
			}
		}
		cris.returnList(ln);
		
		readsOut+=singlesOut+pairsOut;
		basesOut+=singleBasesOut+pairBasesOut;
	}
	
	private void process3_fixInterleaving(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros, final ConcurrentReadOutputStream rosb){

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=ln.list;
		
		final ArrayList<Read> pairs=(ros==null ? null : new ArrayList<Read>(Shared.bufferLen()));
		final ArrayList<Read> singles=(rosb==null ? null : new ArrayList<Read>(Shared.bufferLen()));
		
		Read current=null, prev=null;
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			for(int i=0; i<reads.size(); i++){
				
				current=reads.get(i);
//				if(verbose){System.err.println("Fetched "+current);}
				
				if(prev!=null){
					boolean b=FASTQ.testPairNames(prev, current, allowIdenticalPairNames);
					if(b){
						if(verbose){System.err.println("A");}
						processPair(prev, current, pairs, singles);
						prev=null;
						current=null;
					}else{
						if(verbose){System.err.println("B");}
						processPair(prev, null, null, singles);
						prev=null;
					}
				}
				prev=current;
				current=null;
			}
			
//			if(verbose){System.err.println("X\n"+current+"\n"+prev+"\n");}
			
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
			
			if((ln==null || reads==null || reads.isEmpty()) && prev!=null){ //Process last read
				boolean b=FASTQ.testPairNames(prev, current, allowIdenticalPairNames);
				if(b){
					if(verbose){System.err.println("C");}
					processPair(prev, current, pairs, singles);
					prev=null;
					current=null;
				}else{
					if(verbose){System.err.println("D");}
					processPair(prev, null, null, singles);
					prev=null;
				}
			}
			
			if(rosb!=null){
				if(verbose){System.err.println("Adding "+singles.size()+" to single out.");}
				rosb.add(new ArrayList<Read>(singles), ln.id);
				singles.clear();
			}
			
			if(ros!=null){
				if(verbose){System.err.println("Adding "+pairs.size()+" to pair out.");}
				ros.add(new ArrayList<Read>(pairs), ln.id);
				pairs.clear();
			}
		}
		cris.returnList(ln);
		
		readsOut+=singlesOut+pairsOut;
		basesOut+=singleBasesOut+pairBasesOut;
	}
	
	private void process3_repair(final DualCris cris, final ConcurrentReadOutputStream ros, final ConcurrentReadOutputStream rosb){

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=ln.list;
		
		final ArrayList<Read> pairs=(ros==null ? null : new ArrayList<Read>(Shared.bufferLen()));
		
		boolean foundR1=false, foundR2=false;
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			for(Read r1 : reads){
				Read r2=r1.mate;
				
				if(r1.pairnum()==0){foundR1=true;}
				else{foundR2=true;}
				if(r2!=null){
					if(r2.pairnum()==0){foundR1=true;}
					else{foundR2=true;}
				}
				
				{
					Read pair=repair(r1);
					if(pair!=null && pairs!=null){pairs.add(pair);}
				}
				{
					Read pair=repair(r2);
					if(pair!=null && pairs!=null){pairs.add(pair);}
				}
			}
			
//			if(verbose){System.err.println("X\n"+current+"\n"+prev+"\n");}
			
			cris.returnList(ln.id, foundR1, foundR2);
			foundR1=foundR2=false;
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
			
			if(ros!=null){
				if(verbose){System.err.println("Adding "+pairs.size()+" to pair out.");}
				ros.add(new ArrayList<Read>(pairs), ln.id);
				pairs.clear();
			}
		}
		cris.returnList(ln.id, foundR1, foundR2);
		
		if(!pairMap.isEmpty()){
			final ArrayList<Read> singles=new ArrayList<Read>(pairMap.size());
			for(String key : pairMap.keySet()){
				Read r=pairMap.get(key);
				singles.add(r);
				singlesOut++;
				singleBasesOut+=r.length();
			}
			pairMap.clear();
			if(verbose){System.err.println("Adding "+singles.size()+" to single out.");}
			if(rosb!=null){rosb.add(singles, 0);}
		}
		
		readsOut+=singlesOut+pairsOut;
		basesOut+=singleBasesOut+pairBasesOut;
	}
	
	private void process3_repair(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros, final ConcurrentReadOutputStream rosb){

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=ln.list;
		
		final ArrayList<Read> pairs=(ros==null ? null : new ArrayList<Read>(Shared.bufferLen()));
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			for(Read r1 : reads){
				Read r2=r1.mate;
				
				{
					Read pair=repair(r1);
					if(pair!=null && pairs!=null){pairs.add(pair);}
				}
				{
					Read pair=repair(r2);
					if(pair!=null && pairs!=null){pairs.add(pair);}
				}
			}
			
//			if(verbose){System.err.println("X\n"+current+"\n"+prev+"\n");}
			
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
			
			if(ros!=null && pairs!=null){//pairs!=null is implied
				if(verbose){System.err.println("Adding "+pairs.size()+" to pair out.");}
				ros.add(new ArrayList<Read>(pairs), ln.id);
				pairs.clear();
			}
		}
		cris.returnList(ln);
		
		if(!pairMap.isEmpty()){
			final ArrayList<Read> singles=new ArrayList<Read>(pairMap.size());
			for(String key : pairMap.keySet()){
				Read r=pairMap.get(key);
				singles.add(r);
				singlesOut++;
				singleBasesOut+=r.length();
			}
			pairMap.clear();
			if(verbose){System.err.println("Adding "+singles.size()+" to single out.");}
			if(rosb!=null){rosb.add(singles, 0);}
		}
		
		readsOut+=singlesOut+pairsOut;
		basesOut+=singleBasesOut+pairBasesOut;
	}
	
	
	private int processPair(Read r1, Read r2, ArrayList<Read> pairs, ArrayList<Read> singles){
		int removed=0;
		readsIn++;
		basesIn+=r1.length();
		if(r2!=null){
			readsIn++;
			basesIn+=r2.length();
		}
		
		if(qtrimLeft || qtrimRight){
			if(r1!=null){
				int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, trimE, 1);
				basesTrimmed+=x;
				readsTrimmed+=(x>0 ? 1 : 0);
				if(addcolon){
					String colon=colon1;
					if(!r1.id.contains(colon)){r1.id+=colon;}
				}
			}
			if(r2!=null){
				int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, trimE, 1);
				basesTrimmed+=x;
				readsTrimmed+=(x>0 ? 1 : 0);
				if(addcolon){
					String colon=colon2;
					if(!r2.id.contains(colon)){r2.id+=colon;}
				}
			}
		}
		final int rlen1=(r1==null ? -1 : r1.length());
		final int rlen2=(r2==null ? -1 : r2.length());
		if(verbose){System.err.println("rlen="+rlen1+", rlen2="+rlen2);}
		
		if(rlen1>=minReadLength && rlen2>=minReadLength){
			if(verbose){System.err.println("Sending to pair out:\t"+r1.id+"\t"+r2.id);}
			r1.mate=r2;
			r2.mate=r1;
			r1.setPairnum(0);
			r2.setPairnum(1);
			if(pairs!=null){pairs.add(r1);}
			pairsOut+=2;
			pairBasesOut+=(rlen1+rlen2);
		}else if(rlen1>=minReadLength){
			if(verbose){System.err.println("Sending r1 to single out:\t"+r1.id+"\t"+(r2==null ? "*" : r2.id));}
			r1.mate=null;
			r1.setPairnum(0);
			if(singles!=null){singles.add(r1);}
			singlesOut++;
			singleBasesOut+=rlen1;
			if(r2!=null){removed++;}
		}else if(rlen2>=minReadLength){
			if(verbose){System.err.println("Sending r2 to single out:\t"+(r1==null ? "*" : r1.id)+"\t"+r2.id);}
			r2.mate=null;
			r2.setPairnum(0);
			if(singles!=null){singles.add(r2);}
			singlesOut++;
			singleBasesOut+=rlen2;
			if(r1!=null){removed++;}
		}else{
			if(verbose){System.err.println("Removed both reads:\t"+(r1==null ? "*" : r1.id)+"\t"+(r2==null ? "*" : r2.id));}
			if(r1!=null){removed++;}
			if(r2!=null){removed++;}
		}
		return removed;
	}
	
	
	private Read repair(Read r){
		if(r==null){return null;}
		r.mate=null;
		
		readsIn++;
		basesIn+=r.length();
		final String id=r.id;
		
		final SamLine sl=r.samline;
		if(sl!=null && (!sl.primary() || sl.supplementary())){return null;}
		
		assert(id!=null) : "Read number "+r.numericID+" has no name and thus cannot be re-paired.  To ignore this, run with the -da flag.";
		if(id==null){return null;}
		final int slash=id.indexOf('/');
		String[] split=id.split("\\s+");
		
		if(split.length==1 && slash>0){
			split=new String[] {id.substring(0, slash), id.substring(slash)};
		}
		
		assert(split.length>0);
		String prefix=split[0];
		String suffix=(split.length==1 ? null : split[split.length-1]);
		
		if(sl!=null){
			r.setPairnum(sl.pairnum());
		}else if(suffix!=null){
			if(suffix.startsWith("/1") || suffix.startsWith("1:")){
				r.setPairnum(0);
			}else if(suffix.startsWith("/2") || suffix.startsWith("2:")){
				r.setPairnum(1);
			}else if(id.contains("/1") || id.contains("/2")){
				split=id.split("/");
				prefix=split[0];
				suffix=(split.length==1 ? null : split[split.length-1]);
				
				if(suffix!=null){
					if(suffix.startsWith("1")){
						r.setPairnum(0);
					}else if(suffix.startsWith("2")){
						r.setPairnum(1);
					}
				}else{
					//pairnum cannot be determined
				}
			}else{
				//pairnum cannot be determined
			}
		}else{
			//pairnum cannot be determined
		}
		
		if(addcolon){
			String colon=(r.pairnum()==0 ? colon1 : colon2);
			if(!r.id.contains(colon)){r.id+=colon;}
		}
		
		Read old=pairMap.remove(prefix);
		
//		System.out.println("Processing:\n"+r+"\n"+old+"\n"+readsIn+", "+readsOut+", "+pairsOut);
		
		if(old==null){
			pairMap.put(prefix, r);
			return null;
		}else{
			r.mate=old;
			old.mate=r;
			
			int len=r.length()+old.length();
			pairsOut+=2;
			pairBasesOut+=len;
			
			if(old.pairnum()==1){
				r.setPairnum(0);
				return r;
			}else{
				old.setPairnum(0);
				r.setPairnum(1);
				return old;
			}
		}
	}
	
	
	private String in1=null, in2=null;
	private String out1=null, out2=null;
	private String outsingle=null;
	private long maxReads=-1;
	public boolean errorState=false;
	
	long readsIn=0;
	long basesIn=0;
	long readsOut=0;
	long basesOut=0;
	long pairsOut=0;
	long pairBasesOut=0;
	long singlesOut=0;
	long singleBasesOut=0;
	long readsTrimmed=0;
	long basesTrimmed=0;

	private final LinkedHashMap<String, Read> pairMap;

	private float trimq=6;
	/** Error rate for trimming (derived from trimq) */
	private final float trimE;
	private int minReadLength=20;
	private final boolean qtrimLeft, qtrimRight;
	
	private final boolean fixInterleaving;
	private final boolean allowIdenticalPairNames;
	private final boolean repair;

	private boolean addslash=false;
	private boolean addcolon=false;
	
	private static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	public static boolean showSpeed=true;
	public static boolean verbose=false;
	
	private static final String slash1=" /1";
	private static final String slash2=" /2";
	private static final String colon1=" 1:";
	private static final String colon2=" 2:";
	
}
