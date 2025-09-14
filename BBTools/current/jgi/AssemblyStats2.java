package jgi;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import dna.FastaToChromArrays2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import json.JsonObject;
import shared.KillSwitch;
import shared.MetadataWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.IntList;
import structures.LongList;


/**
 * @author Brian Bushnell
 * @date Sep 16, 2013
 *
 */
public final class AssemblyStats2 {
	
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		AssemblyStats2 x=new AssemblyStats2(args);
		x.process();
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	public AssemblyStats2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, outstream, null, false);
			if(Shared.mainClass==null){Shared.mainClass=this.getClass();}
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_UNPIGZ=ReadWrite.USE_PIGZ=true;
		
		int GCBINS_=200;
		int MINSCAF_=0;
		int gchistdecimals_=-1;
		
		for(int i=0; i<args.length; i++){

			if(true){
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
				}else if(arg.contains("=") && (a.equals("in") || a.equals("ref"))){
					in=b;
				}else if(a.equals("gc") || a.equals("gcout")){
					gc=b;
					if(b==null || "summaryonly".equalsIgnoreCase(b) || "none".equalsIgnoreCase(b)){
						gc=null;
					}
				}else if(a.equals("gchist")){
					gchistFile=b;
					if(b==null || "none".equalsIgnoreCase(b)){
						gchistFile=null;
					}
				}else if(a.equals("gchistdecimals")){
					gchistdecimals_=Integer.parseInt(b);
				}else if(a.equals("gcbins")){
					int x=Integer.parseInt(b);
					if(x>0){GCBINS_=Integer.parseInt(b);}
				}else if(a.equals("shist") || a.equals("scaffoldhist")){
					scaffoldHistFile=b;
					if(b==null || "none".equalsIgnoreCase(b)){
						scaffoldHistFile=null;
					}
				}else if(a.equals("out")){
					out=b;
					if(b==null || "summaryonly".equalsIgnoreCase(b) || "none".equalsIgnoreCase(b)){
						out=null;
					}else if("benchmark".equalsIgnoreCase(b)){
						benchmark=true;
						out=null;
						gc=null;
					}
				}else if(a.equals("benchmark")){
					benchmark=Parse.parseBoolean(b);
					if(benchmark){
						out=null;
						gc=null;
					}
				}else if(a.equals("format")){
					if("json".equalsIgnoreCase(b)){FORMAT=8;}
					else{FORMAT=Integer.parseInt(b);}
					if(FORMAT<0 || FORMAT>8){
						throw new RuntimeException("\nUnknown format: "+FORMAT+"; valid values are 1 through 8.\n");
					}
				}else if(a.equals("json")){
					if(Parse.parseBoolean(b)){FORMAT=8;}
					else{FORMAT=1;}
				}else if(a.equals("gcformat")){
					GCFORMAT=Integer.parseInt(b);
					if(GCFORMAT<0 || GCFORMAT>5){
						throw new RuntimeException("\nUnknown gcformat: "+GCFORMAT+"; valid values are 0 through 5.\n");
					}
				}else if(a.equals("gcfraction")){
					GC_FRACTION=Parse.parseBoolean(b);
				}else if(a.equals("cutoff")){
					cutoff=Parse.parseKMG(b);
				}else if(a.equals("k") || a.equals("bbmapkmer")){
					bbmapkmer=Integer.parseInt(b);
				}else if(a.equals("printl90") || a.equals("l90") || a.equals("printn90") || a.equals("n90")){
					printL90=Parse.parseBoolean(b);
				}else if(a.equals("printextended") || a.equals("extended")){
					printExtended=Parse.parseBoolean(b);
				}else if(a.equals("overwrite") || a.equals("ow")){
					overwrite=Parse.parseBoolean(b);
				}else if(a.equals("n_")){
					N_UNDERSCORE=Parse.parseBoolean(b);
				}else if(a.equals("header") || a.equals("useheader")){
					useheader=Parse.parseBoolean(b);
				}else if(a.equals("addfilename") || a.equals("addname")){
					addfilename=Parse.parseBoolean(b);
				}else if(a.equals("minscaf") || a.equals("mincontig") || a.equals("minlen") || a.equals("min")){
					MINSCAF_=Parse.parseIntKMG(b);
				}else if(a.equals("showspeed") || a.equals("ss")){
					showspeed=Parse.parseBoolean(b);
				}else if(a.equals("printheadersize") || a.equals("phs")){
					printheadersize=Parse.parseBoolean(b);
				}else if(a.equals("skipduplicatelines") || a.equals("sdl")){
					skipDuplicateLines=Parse.parseBoolean(b);
				}else if(a.equals("printduplicatelines") || a.equals("pdl")){
					skipDuplicateLines=!Parse.parseBoolean(b);
				}else if(a.equals("showbbmap")){
					if(!Parse.parseBoolean(b)){bbmapkmer=0;}
				}else if(a.equals("contigbreak") || (arg.contains("=") && (a.equals("n")))){
					maxNs=Integer.parseInt(b);
				}
				
				else if(a.equals("logsumoffset") || a.equals("logoffset")){
					logSumOffset=Integer.parseInt(b);
				}else if(a.equals("logsumbase") || a.equals("logbase")){
					logSumBase=Double.parseDouble(b);
				}else if(a.equals("squarelog") || a.equals("logsquare") || a.equals("logsquared")){
					squareLog=Parse.parseBoolean(b);
				}else if(a.equals("logpower")){
					logPower=Double.parseDouble(b);
				}else if(a.equals("powsum") || a.equals("powersum") || a.equals("powersumpower") || a.equals("powsumpower")){
					powSumPower=Double.parseDouble(b);
				}
				
				else if(a.equalsIgnoreCase("assemblyScoreMinLen")){
					assemblyScoreMinLen=Parse.parseKMG(b);
				}else if(a.equalsIgnoreCase("assemblyScoreMaxLen")){
					assemblyScoreMaxLen=Parse.parseKMG(b);
				}else if(a.equalsIgnoreCase("alignedFraction") || a.equalsIgnoreCase("aligned") || a.equalsIgnoreCase("mapped")){
					alignedFraction=Double.parseDouble(b);
					if(alignedFraction>1){alignedFraction/=100;}
				}else if(a.equals("score") || a.equals("assemblyscore")){
					printAssemblyScore=Parse.parseBoolean(b);
				}
				
				else if(i>0 && (a.equals("n") || a.equals("-n")) && b!=null){
					maxNs=Integer.parseInt(b);
				}else if(in==null && i==0 && !arg.contains("=")){
					in=arg;
				}else{
					throw new RuntimeException("Unknown parameter "+arg);
				}
			}
		}
		
		if(printExtended){
			printL90=true;
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		minScaffold=MINSCAF_;
		gcbins=GCBINS_;
		
		if(gchistdecimals_<1){
			gchistdecimals_=3;
			if(gcbins==2 || gcbins==5 || gcbins==10){
				gchistdecimals_=1;
			}else if(gcbins==20 || gcbins==25 || gcbins==50 || gcbins==100){
				gchistdecimals_=2;
			}
			if(gcbins>1000 && gchistdecimals_<4){gchistdecimals_=4;}
			if(gcbins>10000 && gchistdecimals_<5){gchistdecimals_=5;}
		}
		gchistDecimals1=gchistdecimals_;
		
		if(maxNs<0){maxNs=10;}

		if(out==null || out.equalsIgnoreCase("stdout") || out.equalsIgnoreCase("standardout")){out=null;}
		
		clist=new LongList((int)Tools.min(1<<15, cutoff+1)); //Number of contigs of length x
		slist=new LongList((int)Tools.min(1<<15, cutoff+1)); //Number of scaffolds of length x
		sclist1=new LongList((int)Tools.min(1<<15, cutoff+1)); //Sum of contigs per scaffold of length x
		sclist2=new LongList((int)Tools.min(1<<15, cutoff+1)); //Sum of contig lengths per scaffold of length x
		
		llist=new LongList(64); //List of contig lengths for contigs at least cutoff in length
		tlist=new ArrayList<Triple>(64); //List of scaf len, contigs, contig sum for scaffolds at least cutoff in length
		
		gcbins2=(gcbins>=1000 ? gcbins : gcbins*10);
		
		gchistArray=new long[gcbins2];
		gchist_by_base=new long[gcbins2];
		
	}
	
	/*--------------------------------------------------------------*/
	
	public void process(){
		Timer t=new Timer();
		
		InputStream is=null;
		{
			if(in==null){throw new RuntimeException("No input file.");}
			if(in.equalsIgnoreCase("stdin") || in.equalsIgnoreCase("standardin")){
				is=System.in;
			}else{
				File f=new File(in);
				if((!f.exists() || f.isDirectory()) && !in.toLowerCase().startsWith("stdin")){
					throw new RuntimeException("Input file does not appear to be valid: "+in);
				}
			}
		}
		
		long[] counts=null;
		long sum=0;
		
		boolean fastqMode=false;
		if(is!=System.in){
			FileFormat ff=null;
			try {
				ff=FileFormat.testInput(in, FileFormat.FA, null, false, true, true);
			} catch (Throwable e) {
				e.printStackTrace();
			}
			if(ff!=null){
				fastqMode=ff.fastq();
			}
			assert(!ff.samOrBam()) : "AssemblyStats does not support sam/bam files, only fasta and fastq.";
//			assert(ff==null || (!ff.fastq())) : "AssemblyStats only supports fasta files.  To override this message, use the -da flag.";
		}
		
		if(is==null){is=ReadWrite.getInputStream(in, false, true);}
		try {
			if(benchmark){sum=bench(is);}
			else{
				if(fastqMode){
					counts=countFastq(is, gc);
				}else{
					counts=countFasta(is, gc);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			if(is!=System.in){is.close();}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		if(tlist!=null && tlist.size()>0){Shared.sort(tlist);}
		if(llist!=null && llist.size>0){Shared.sort(llist.array, 0, llist.size);}
		
		t.stop();
		
		if(benchmark){
			printBenchResults(t, counts, sum, in);
		}else{
//			outstream.println("\nclist="+clist+"\nslist="+slist+"\nsclist1="+sclist1+"\nsclist2="+sclist2+"\nllist="+llist+"\ntlist="+tlist+"\n"); //***
//			outstream.println("\nclist.size="+clist.size+"\nslist.size="+slist.size+"\nsclist1.size="+sclist1.size+"\nsclist2.size="+sclist2.size+"\nllist.size="+llist.size+"\ntlist.size()="+tlist.size()+"\n"); //***
			printResults(t, counts, sum, minScaffold, gc_std, in, clist, slist, sclist1, sclist2, llist, tlist, out);
//			outstream.println("Printed results to "+out);
			writeHistFile(scaffoldHistFile, slist, tlist, false);
			
			if(gchistFile!=null){printGCHist(gchistFile);}
		}
	}
	
	/*--------------------------------------------------------------*/
	
	public long[] countFasta(final InputStream is, String gcout) throws IOException{
		
		long limsum=0;
		long headerlen=0;
		final byte[] buf=new byte[32768];
		final TextStreamWriter tswgc=(gcout==null ? null : new TextStreamWriter(gcout, overwrite, append, false));
		if(tswgc!=null){
			tswgc.start();
			if(GCFORMAT==0 || GCFORMAT==1){
				tswgc.println("#Name\tLength\tA\tC\tG\tT\tN\tIUPAC\tOther\tGC");
			}else if(GCFORMAT==2){
				tswgc.println("#Name\tGC");
			}else if(GCFORMAT==3){
				tswgc.println("#Name\tLength\tA\tC\tG\tT\tN\tIUPAC\tOther\tGC");
			}else if(GCFORMAT==4){
				tswgc.println("#Name\tLength\tGC");
			}else if(GCFORMAT==5){
				tswgc.println("#Name\tLength\tGC\tlogsum\tpowersum");
			}else{
				throw new RuntimeException("Unknown format.");
			}
		}
//		assert(false) : GCFORMAT+", "+out+", "+tswgc;
		final long[] counts=KillSwitch.allocLong1D(8);
		final long[] overall=KillSwitch.allocLong1D(8);
		final StringBuilder hdr=(gcout==null ? null : new StringBuilder());
		boolean hdmode=false;
		
		int i=0;
		int lim=is.read(buf);
		limsum+=lim;
		
		int contigs=0;
		int contiglen=0;
//		int contiglensum=0;
		int scaffoldlen=0;
		int ns=0;
		
		final IntList currentContigs=new IntList(10000);
		
		while(lim>0){
			if(hdmode){//Scan to end of header.
				if(hdr==null){//Header is not being stored
					while(i<lim){
						final byte c=buf[i];
						i++;
						if(c<=slashr){
							hdmode=false;
							contiglen=0;
//							contiglensum=0;
							scaffoldlen=0;
							ns=0;
							contigs=0;
							break;
						}
						headerlen++;
					}
				}else{//Store the header
					while(i<lim){
						final byte c=buf[i];
						i++;
						if(c<=slashr){
							hdmode=false;
							contiglen=0;
//							contiglensum=0;
							scaffoldlen=0;
							ns=0;
							contigs=0;
							break;
						}
						hdr.append((char)c);
					}
				}
			}
			
			if(!hdmode){//Scan bases
				while(i<lim){
					final byte c=buf[i];
					assert(c>=0) : i+(hdr==null ? "Non-ascii input." : "\nError with header: "+hdr);
					final byte cnum=charToNum[c];
					i++;

					if(c==carrot){//Start of a new header
//						assert(false) : scaffoldlen;
						hdmode=true;
						if(scaffoldlen>0){//if the scaffold was not blank
							
							{//NEW
								if(contiglen>0 || contigs==0){
									currentContigs.set(contigs, contiglen);
									contigs++;
//									System.out.println("For header "+hdr+": added contig.  len="+contiglen+", contigs="+contigs);
//									contiglensum+=contiglen;
								}
							}
							
//							assert(false);
							if(scaffoldlen>=minScaffold){

								int contiglensum=0;
								{//NEW
//									System.out.println("Dumping "+contigs+" contigs.");
									for(int j=0; j<contigs; j++){
										final int cl=currentContigs.get(j);
										if(cl>0 || contigs==0){
											contiglensum+=cl;
											if(cl<cutoff){
												clist.increment(cl, 1);
											}else{
												llist.add(cl);
											}
										}
									}
								}
								
								if(scaffoldlen<cutoff){
									slist.increment(scaffoldlen, 1);
									sclist1.increment(scaffoldlen, contigs);
									sclist2.increment(scaffoldlen, contiglensum);
								}else{
									tlist.add(new Triple(scaffoldlen, contigs, contiglensum));
								}


								if(hdr!=null){
									if(tswgc!=null){tswgc.print(toString2(hdr, counts));}//if added for compiler benefit
									headerlen+=hdr.length();
									hdr.setLength(0);
								}
								{
									long gc=counts[1]+counts[2];
									long acgt=gc+counts[0]+counts[3];
									if(acgt>0){
										int index=Tools.min((int)((gc*gcbins2)/acgt),gcbins2-1);
										gchistArray[index]++;
										gchist_by_base[index]+=scaffoldlen;
//										assert(false);
									}
								}
								for(int j=0; j<counts.length; j++){
									overall[j]+=counts[j];
									counts[j]=0;
								}
							}else{
								Arrays.fill(counts, 0);
								if(hdr!=null){hdr.setLength(0);}
							}
						}

						break;
					}
					
					if(c>slashr){
						counts[cnum]++;
						scaffoldlen++;
						
//						if(c!=noref && c!=noref2){
						if(cnum!=5){
							ns=0;
							contiglen++;
						}else{
							ns++;
							if(ns==maxNs && contiglen>0){
//								if(contiglen<cutoff){
//									clist.increment(contiglen, 1);
//								}else{
//									llist.add(contiglen);
//								}
////								clist.increment(contiglen, 1);
//								contiglensum+=contiglen;
//								contiglen=0;
//								contigs++;

								{//NEW
									currentContigs.set(contigs, contiglen);
									contiglen=0;
									contigs++;
								}
							}
						}
					}
				}
			}
			if(i>=lim){
				i=0;
				lim=is.read(buf);
				limsum+=lim;
			}
		}

		if(scaffoldlen>0){
			
//			if(contiglen>0 || contigs==0){
//				contigs++;
//				contiglensum+=contiglen;
//				if(contiglen<cutoff){
//					clist.increment(contiglen, 1);
//				}else{
//					llist.add(contiglen);
//				}
//			}
			
			{//NEW
				if(contiglen>0 || contigs==0){
					currentContigs.set(contigs, contiglen);
					contigs++;
//					contiglensum+=contiglen;
				}
			}
			
			if(scaffoldlen>=minScaffold){

				int contiglensum=0;
				{//NEW
//					System.out.println("Dumping "+contigs+" contigs.");
					for(int j=0; j<contigs; j++){
						final int cl=currentContigs.get(j);
						if(cl>0 || contigs==0){
							contiglensum+=cl;
							if(cl<cutoff){
								clist.increment(cl, 1);
							}else{
								llist.add(cl);
							}
						}
					}
				}
				
				if(scaffoldlen<cutoff){
					slist.increment(scaffoldlen, 1);
					sclist1.increment(scaffoldlen, contigs);
					sclist2.increment(scaffoldlen, contiglensum);
				}else{
					tlist.add(new Triple(scaffoldlen, contigs, contiglensum));
				}


//				slist.increment(scaffoldlen, 1);
//				if(contiglen>0 || contigs==0){
//					contigs++;
//					contiglensum+=contiglen;
//					clist.increment(contiglen, 1);
//				}
//				sclist1.increment(scaffoldlen, contigs);
//				sclist2.increment(scaffoldlen, contiglensum);

				if(hdr!=null){
					if(tswgc!=null){tswgc.print(toString2(hdr, counts));}//if added for compiler benefit
					hdr.setLength(0);
				}

				{
					long gc=counts[1]+counts[2];
					long acgt=gc+counts[0]+counts[3];
					if(acgt>0){
						int index=Tools.min((int)((gc*gcbins2)/acgt),gcbins2-1);
						gchistArray[index]++;
						gchist_by_base[index]+=scaffoldlen;
					}
				}
				for(int j=0; j<counts.length; j++){
					overall[j]+=counts[j];
					counts[j]=0;
				}
			}
		}
		
//		outstream.println("clist="+clist+"\nslist="+slist+"\nsclist1="+sclist1+"\nsclist2="+sclist2+"\nllist="+llist+"\ntlist="+tlist); //***
		
		
		if(tswgc!=null){
			if(tswgc.fname.equalsIgnoreCase("stdout") || tswgc.fname.startsWith("stdout.")){
				if(FORMAT>0 && (out==null || out.equalsIgnoreCase("stdout") || out.startsWith("stdout."))){
					tswgc.print("\n");
				}
			}
			tswgc.poison();
			tswgc.waitForFinish();
		}
		LIMSUM=limsum;
		HEADERLENSUM=headerlen;
		

		gc_std=Tools.standardDeviationHistogram(gchistArray)/gcbins2;
		gchistArray_downsampled=Tools.downsample(gchistArray, gcbins);

		gc_bb_std=Tools.standardDeviationHistogram(gchist_by_base)/gcbins2;
		gchist_by_base_downsampled=Tools.downsample(gchist_by_base, gcbins);
		
		return overall;
	}
	
	public long[] countFastq(final InputStream is, String gcout) throws IOException{
		
		long limsum=0;
		long headerlen=0;
		final byte[] buf=new byte[32768];
		final TextStreamWriter tswgc=(gcout==null ? null : new TextStreamWriter(gcout, overwrite, false, false));
		if(tswgc!=null){
			tswgc.start();
			if(GCFORMAT==0 || GCFORMAT==1){
				tswgc.println("#Name\tLength\tA\tC\tG\tT\tN\tIUPAC\tOther\tGC");
			}else if(GCFORMAT==2){
				tswgc.println("#Name\tGC");
			}else if(GCFORMAT==3){
				tswgc.println("#Name\tLength\tA\tC\tG\tT\tN\tIUPAC\tOther\tGC");
			}else if(GCFORMAT==4){
				tswgc.println("#Name\tLength\tGC");
			}else if(GCFORMAT==5){
				tswgc.println("#Name\tLength\tGC\tlogsum\tpowersum");
			}else{
				throw new RuntimeException("Unknown format.");
			}
		}
//		assert(false) : GCFORMAT+", "+out+", "+tswgc;
		final long[] counts=KillSwitch.allocLong1D(8);
		final long[] overall=KillSwitch.allocLong1D(8);
		final StringBuilder hdr=(gcout==null ? null : new StringBuilder());
		int line=0;
		
		int i=0;
		int lim=is.read(buf);
		limsum+=lim;
		
		int contigs=0;
		int contiglen=0;
//		int contiglensum=0;
		int scaffoldlen=0;
		int ns=0;
		
		final IntList currentContigs=new IntList(10000);
		
		while(lim>0){
			if(line==0){//Scan to end of header.
//				outstream.println("1");
				if(hdr==null){//Header is not being stored
					while(i<lim){
						final byte c=buf[i];
						i++;
						if(c<=slashr){
							line++;
//							outstream.println("1.1");
							contiglen=0;
//							contiglensum=0;
							scaffoldlen=0;
							ns=0;
							contigs=0;
							break;
						}
						headerlen++;
					}
				}else{//Store the header
					while(i<lim){
						final byte c=buf[i];
						i++;
						if(c<=slashr){
							line++;
//							outstream.println("1.2");
							contiglen=0;
//							contiglensum=0;
							scaffoldlen=0;
							ns=0;
							contigs=0;
							break;
						}
						hdr.append((char)c);
					}
				}
			}
			
			if(line==1){//Scan bases
//				outstream.println("2");
				while(i<lim){
					final byte c=buf[i];
					final byte cnum=charToNum[c];
					i++;

					if(c<=slashr){//Finish the contig
//						assert(false) : scaffoldlen;
						line=(line+1)&3;
//						outstream.println("2.1");
						if(scaffoldlen>0){//if the scaffold was not blank
							
							{//NEW
								if(contiglen>0 || contigs==0){
									currentContigs.set(contigs, contiglen);
									contigs++;
//									System.out.println("For header "+hdr+": added contig.  len="+contiglen+", contigs="+contigs);
//									contiglensum+=contiglen;
								}
							}
							
//							assert(false);
							if(scaffoldlen>=minScaffold){

								int contiglensum=0;
								{//NEW
//									System.out.println("Dumping "+contigs+" contigs.");
									for(int j=0; j<contigs; j++){
										final int cl=currentContigs.get(j);
										if(cl>0 || contigs==0){
											contiglensum+=cl;
											if(cl<cutoff){
												clist.increment(cl, 1);
											}else{
												llist.add(cl);
											}
										}
									}
								}
								
								if(scaffoldlen<cutoff){
									slist.increment(scaffoldlen, 1);
									sclist1.increment(scaffoldlen, contigs);
									sclist2.increment(scaffoldlen, contiglensum);
								}else{
									tlist.add(new Triple(scaffoldlen, contigs, contiglensum));
								}


								if(hdr!=null){
									if(tswgc!=null){tswgc.print(toString2(hdr, counts));}//if added for compiler benefit
									headerlen+=hdr.length();
									hdr.setLength(0);
								}
								{
									long gc=counts[1]+counts[2];
									long acgt=gc+counts[0]+counts[3];
									if(acgt>0){
										int index=Tools.min((int)((gc*gcbins2)/acgt),gcbins2-1);
										gchistArray[index]++;
										gchist_by_base[index]+=scaffoldlen;
//										assert(false);
									}
								}
								for(int j=0; j<counts.length; j++){
									overall[j]+=counts[j];
									counts[j]=0;
								}
							}else{
								Arrays.fill(counts, 0);
								if(hdr!=null){hdr.setLength(0);}
							}
						}

						break;
					}
					
					if(c>slashr){
						counts[cnum]++;
						scaffoldlen++;
						
//						if(c!=noref && c!=noref2){
						if(cnum!=5){
							ns=0;
							contiglen++;
						}else{
							ns++;
							if(ns==maxNs && contiglen>0){
//								if(contiglen<cutoff){
//									clist.increment(contiglen, 1);
//								}else{
//									llist.add(contiglen);
//								}
////								clist.increment(contiglen, 1);
//								contiglensum+=contiglen;
//								contiglen=0;
//								contigs++;

								{//NEW
									currentContigs.set(contigs, contiglen);
									contiglen=0;
									contigs++;
								}
							}
						}
					}
				}
			}

//			outstream.println("3");
			if(i>=lim){
				i=0;
				lim=is.read(buf);
				limsum+=lim;
			}

//			outstream.println("4");
			assert(line>1 || lim<=i || i==0) : line+", "+i+", "+lim;
			while(i<lim && line>1){
				final byte c=buf[i];
				i++;
				if(c<=slashr){
					line=(line+1)&3;
				}
			}
		}
		
		if(tswgc!=null){
			if(tswgc.fname.equalsIgnoreCase("stdout") || tswgc.fname.startsWith("stdout.")){
				if(FORMAT>0 && (out==null || out.equalsIgnoreCase("stdout") || out.startsWith("stdout."))){
					tswgc.print("\n");
				}
			}
			tswgc.poison();
			tswgc.waitForFinish();
		}
		LIMSUM=limsum;
		HEADERLENSUM=headerlen;
		

		gc_std=Tools.standardDeviationHistogram(gchistArray)/gcbins2;
		gchistArray_downsampled=Tools.downsample(gchistArray, gcbins);

		gc_bb_std=Tools.standardDeviationHistogram(gchist_by_base)/gcbins2;
		gchist_by_base_downsampled=Tools.downsample(gchist_by_base, gcbins);
		
		return overall;
	}
	
	/*--------------------------------------------------------------*/
	
	private void printGCHist(String gchistFile){
		if(!Tools.canWrite(gchistFile, overwrite)){
			outstream.println("Can't write gc histogram because file exists and overwrite="+overwrite);
			assert(false);
		}else{
			long gchistFilesum=shared.Vector.sum(gchistArray_downsampled);
			long gchistFilesumbb=shared.Vector.sum(gchist_by_base_downsampled);
			double invsum=(gchistFilesum==0 ? 0 : 1.0/gchistFilesum);
			double invsumbb=(gchistFilesum==0 ? 0 : 1.0/gchistFilesumbb);
			double invbins=1.0/(gcbins==0 ? 1 : gcbins);
//			assert(false) : Arrays.toString(gchistArray);
			StringBuilder sb=new StringBuilder();
			sb.append(Tools.format("#GC\tscaffolds\tfraction\tlength\tlen_fraction\n"));
			for(int i=0; i<gcbins; i++){
				sb.append(Tools.format("%."+gchistDecimals1+"f\t%d\t%.5f\t%d\t%.5f\n",
						i*invbins, gchistArray_downsampled[i], gchistArray_downsampled[i]*invsum, gchist_by_base_downsampled[i], gchist_by_base_downsampled[i]*invsumbb));
			}
			if(gchistFile.equalsIgnoreCase("stdout")){
				System.out.println(sb);
			}else{
				ReadWrite.writeString(sb, gchistFile);
			}
		}
	}
	
	
	public static void printBenchResults(Timer t, long[] counts, long sum, String in){
		outstream.println("Time: \t"+t);
		long bytes=new File(in).length();
		if(bytes<1){bytes=LIMSUM;}
		double mbps1=bytes*1000d/t.elapsed;
		double mbps2=sum*1000d/t.elapsed;
		outstream.println(Tools.format("Raw Speed:         \t%.2f MBytes/s",mbps1));
		outstream.println(Tools.format("Uncompressed Speed:\t%.2f MBytes/s",mbps2));
	}
	
	public static double calcLogSumContigs(LongList clist, LongList llist, int cutoff, double base){
		return calcLogSumCounts(clist, cutoff, base)+calcLogSumLengths(llist, cutoff, base);
	}
	
	public static double calcLogSumScaffolds(LongList slist, ArrayList<Triple> tlist, int cutoff, double base){
		return calcLogSumCounts(slist, cutoff, base)+calcLogSumTriples(tlist, cutoff, base);
	}
	
	public static double calcPowerSumContigs(LongList clist, LongList llist, int cutoff, double power){
		return calcPowerSumCounts(clist, cutoff, power)+calcPowerSumLengths(llist, cutoff, power);
	}
	
	public static double calcPowerSumScaffolds(LongList slist, ArrayList<Triple> tlist, int cutoff, double power){
		return calcPowerSumCounts(slist, cutoff, power)+calcPowerSumTriples(tlist, cutoff, power);
	}
	
	public static double calcAssemblyScore(LongList slist, ArrayList<Triple> tlist, long min, long max, double aligned){
		long lengthSum=0;
		double scoreSum=0;
		for(int i=(int)min; i<slist.size(); i++){
			long count=slist.get(i);
			lengthSum+=(i*count);
		}
		for(Triple t : tlist){
			lengthSum+=t.length;
		}
		
		double mult=1.0/(1000*lengthSum);
		for(int i=(int)min; i<slist.size(); i++){
			long count=slist.get(i);
			if(count>0){
				scoreSum+=count*(i*Tools.min(max, i))*mult;
			}
		}
		for(Triple t : tlist){
			scoreSum+=(t.length*Tools.min(max, t.length))*mult;
		}
		
		return scoreSum+aligned*50;
	}
	
	public static double calcLogSumCounts(LongList counts, int cutoff, double base){
		final double mult=1/Math.log(base);
		
		double sum=0;
		for(int i=cutoff, max=counts.size(); i<max; i++){
			long count=counts.get(i);
			if(count>0){
				double log=mult*Math.log(i);
				if(squareLog){log=log*log*0.05;}
				else if(logPower!=1){log=Math.pow(log*0.12, logPower);}
				double incr=log*count*i;
				sum+=incr;
			}
		}
		return sum;
	}
	
	private static double logsum(long len){
		if(len<logSumOffset){return 0;}
		final double mult=1/Math.log(logSumBase);
		double log=mult*Math.log(len);
		if(squareLog){log=log*log*0.05;}
		else if(logPower!=1){log=Math.pow(log*0.12, logPower);}
		double incr=log*len;
		return incr/logSumOffset;
	}
	
	public static double calcLogSumLengths(LongList lengths, int cutoff, double base){
		final double mult=1/Math.log(base);
		
		double sum=0;
		for(int i=0, max=lengths.size(); i<max; i++){
			long length=lengths.get(i);
			if(length>=cutoff){
				double log=mult*Math.log(length);
				if(squareLog){log=log*log*0.05;}
				else if(logPower!=1){log=Math.pow(log*0.12, logPower);}
				double incr=log*length;
				sum+=incr;
			}
		}
		return sum;
	}
	
	public static double calcLogSumTriples(ArrayList<Triple> triples, int cutoff, double base){
		final double mult=1/Math.log(base);
		
		double sum=0;
		for(int i=0, max=triples.size(); i<max; i++){
			long length=triples.get(i).length;
			if(length>=cutoff){
				double log=mult*Math.log(length);
				if(squareLog){log=log*log*0.05;}
				else if(logPower!=1){log=Math.pow(log*0.12, logPower);}
				double incr=log*length;
				sum+=incr;
			}
		}
		return sum;
	}
	
	public static double calcPowerSumCounts(LongList counts, int cutoff, double power){
		final double mult=1/Math.pow(1000, power);
		
		double sum=0;
		for(int i=cutoff, max=counts.size(); i<max; i++){
			long count=counts.get(i);
			if(count>0){
				double pow=mult*Math.pow(i, power);
				double incr=pow*count*i;
				sum+=incr;
			}
		}
		return sum;
	}
	
	private static double powersum(long len){
		if(len<logSumOffset){return 0;}
		final double mult=1/Math.pow(1000, powSumPower);
		double pow=mult*Math.pow(len, powSumPower);
		double incr=pow*len;
		return incr/logSumOffset;
	}
	
	public static double calcPowerSumLengths(LongList lengths, int cutoff, double power){
		final double mult=1/Math.pow(1000, power);
		
		double sum=0;
		for(int i=0, max=lengths.size(); i<max; i++){
			long length=lengths.get(i);
			if(length>=cutoff){
				double pow=mult*Math.pow(length, power);
				double incr=pow*length;
				sum+=incr;
			}
		}
		return sum;
	}
	
	public static double calcPowerSumTriples(ArrayList<Triple> triples, int cutoff, double power){
		final double mult=1/Math.pow(1000, power);
		
		double sum=0;
		for(int i=0, max=triples.size(); i<max; i++){
			long length=triples.get(i).length;
			if(length>=cutoff){
				double pow=mult*Math.pow(length, power);
				double incr=pow*length;
				sum+=incr;
			}
		}
		return sum;
	}
	
	public static void printResults(Timer t, long[] counts, long sum, long minScaffold, double gc_std, String in, LongList clist, LongList slist, LongList sclist1, LongList sclist2,
			LongList llist, ArrayList<Triple> tlist, String out){
		
		String name=in;
		if(in!=null && !in.toLowerCase().startsWith("stdin")){
			try {
				File f=new File(in);
				name=f.getCanonicalPath();
			} catch (IOException e) {}
		}
		
		long contigs=0;
		long scaffolds=0;
		long contiglen=0;
		long scaflen=0;
		long contigs1;
		long contiglen2;
		long maxScaf=0, maxContig=0;
		long[] carray=clist.array;
		long[] sarray=slist.array;
		long[] scarray1=sclist1.array;
		long[] scarray2=sclist2.array;

		long[] larray=llist.array;
		
		StringBuilder sb=new StringBuilder(), sb2=new StringBuilder();
		
		for(int i=0; i<carray.length; i++){
			long x=carray[i];
			if(x>0){
				contigs+=x;
				contiglen+=(x*i);
				maxContig=i;
			}
		}
		
		for(int i=0; i<sarray.length; i++){
			long x=sarray[i];
			if(x>0){
				scaffolds+=x;
				scaflen+=(x*i);
				maxScaf=i;
			}
		}
		
		contigs+=llist.size;
		for(int i=0; i<llist.size; i++){
			long x=larray[i];
			assert(x>0);
			contiglen+=x;
			maxContig=Tools.max(maxContig, x);
		}
		
		scaffolds+=tlist.size();
		for(Triple tp : tlist){
			scaflen+=tp.length;
			maxScaf=Tools.max(maxScaf, tp.length);
		}
		
		if(FORMAT<3){
			sb.append("Main genome scaffold total:         \t"+scaffolds+"\n");
			sb.append("Main genome contig total:           \t"+contigs+"\n");
		}else if(FORMAT==7){
			sb.append("Main genome contig total:           \t"+contigs+"\n");
		}
		
		if(FORMAT==0){
		}else if(FORMAT==1){
			sb.append("Main genome scaffold sequence total:\t"+Tools.format("%.3f Mb",scaflen/1000000f)+"\n");
			sb.append("Main genome contig sequence total:  \t"+Tools.format("%.3f Mb  \t%.3f%% gap",contiglen/1000000f,(scaflen-contiglen)*100f/scaflen)+"\n");
		}else if(FORMAT==2){
			sb.append("Main genome scaffold sequence total:\t"+scaflen+"\n");
			sb.append("Main genome contig sequence total:  \t"+Tools.format("%d  \t%.3f%% gap",contiglen,(scaflen-contiglen)*100f/scaflen)+"\n");
		}else if(FORMAT==3 || FORMAT==6){
			
		}else if(FORMAT==4){
			
		}else if(FORMAT==5){
			
		}else if(FORMAT==8){
			
		}else if(FORMAT==7){
			sb.append("Main genome contig sequence total:  \t"+Tools.format("%.3f Mb",contiglen/1000000f)+"\n");
		}else{throw new RuntimeException("Unknown format");}
		
		if(FORMAT<3){
			sb2.append("\n");
			sb2.append("Minimum \tNumber        \tNumber        \tTotal         \tTotal         \tScaffold\n");
			sb2.append("Scaffold\tof            \tof            \tScaffold      \tContig        \tContig  \n");
			sb2.append("Length  \tScaffolds     \tContigs       \tLength        \tLength        \tCoverage\n");
			sb2.append("--------\t--------------\t--------------\t--------------\t--------------\t--------\n");
		}else if(FORMAT==7){
			sb2.append("\n");
			sb2.append("Minimum \tNumber        \tTotal         \n");
			sb2.append("Contig  \tof            \tContig        \n");
			sb2.append("Length  \tContigs       \tLength        \n");
			sb2.append("--------\t--------------\t--------------\n");
		}
		
		final int[] lims;
		
		if(FORMAT==7){
			int minScaf=-1;
			for(int i=0; i<carray.length && minScaf<0; i++){
				long x=sarray[i];
				if(x>0){minScaf=i;}
			}
			if(minScaf<0 && !tlist.isEmpty()){
				minScaf=(int)Tools.min(tlist.get(0).length, 250000000);
			}
			if(minScaf<1){minScaf=1;}
			int[] temp=new int[] {0, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000,
					250000, 500000, 1000000, 2500000, 5000000, 10000000, 25000000, 50000000, 100000000, 250000000};
			ArrayList<Integer> tempList=new ArrayList<Integer>();
			tempList.add(0);
			for(int i=1; i<temp.length; i++){
				int x=temp[i];
				int prev=temp[i-1];
				if(x>=minScaf || i==temp.length-1){
					if(prev<minScaf){
						tempList.add(minScaf);
					}else{
						tempList.add(x);
					}
				}
			}
			lims=new int[tempList.size()];
			for(int i=0; i<lims.length; i++){
				lims[i]=tempList.get(i);
			}
		}else{
			lims=new int[] {0, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000,
					250000, 500000, 1000000, 2500000, 5000000, 10000000, 25000000, 50000000, 100000000, 250000000};
		}
		
		int lidx=0;
		int next=0;
		long csum=contigs;
		long ssum=scaffolds;
		long clen=contiglen;
		long slen=scaflen;
		
		long ln50=-1;
		long ll50=-1;
		long cn50=-1;
		long cl50=-1;

		long ln90=-1;
		long ll90=-1;
		long cn90=-1;
		long cl90=-1;

		long s50=slen/2;
		long c50=clen/2;
		
		long s90=(int)(slen*0.9f);
		long c90=(int)(clen*0.9f);

		final int numOverCutoff=(FORMAT==7 ? 1000 : 50000);
		final int numOverCutoffPlusOne=numOverCutoff+1;
		long numOver50=0;
		long basesOver50=0;
		float fractionOver50=0;
		
		
		//Disable printing of 50~500 when not needed
//		{
//			boolean b=true;
//			for(int i=0; i<500 && i<sarray.length && b; i++){
//				b=(sarray[i]==0);
//			}
//			if(b){
//				lidx=Arrays.binarySearch(lims, 1000)-1;
//				lims[lidx]=0;
//			}
//		}
		if(skipDuplicateLines){
			int maxZero=-1;
			for(int i=0; i<sarray.length && sarray[i]==0; i++){
				maxZero=i;
			}
			if(maxZero>=50){
				for(int i=1; i<lims.length; i++){
					if(lims[i]<=maxZero){
						lidx=i-1;
						lims[lidx]=0;
					}else{break;}
				}
			}
		}
		
		final int lim=slist.size;
		assert(lim<=sarray.length);
		
		{
			//These two loops generate the scaffold table
			long prevSsum=-1, prevCsum=-1, prevSlen=-1, prevClen=-1;
			String prevLine=null;
			for(int i=0; i<lim && slen>0; i++){
				//			System.out.println("\n<A1>\ti="+i+", lidx="+lidx+", lims.length="+lims.length+", next="+next+", sarray.length="+sarray.length+", slen="+slen);

				if(i==next){
					prevLine=formatX(next, ssum, csum, slen, clen);
//					if(!skipDuplicateLines || ssum!=prevSsum || csum!=prevCsum || slen!=prevSlen || clen!=prevClen){
					if(prevLine!=null && (next==0 || next>=minScaffold)){
						sb2.append(prevLine);
						sb2.append('\n');
						prevLine=null;
					}
//					}
					prevSsum=ssum; prevCsum=csum; prevSlen=slen; prevClen=clen;
					
					lidx++;
					if(lidx<lims.length){next=lims[lidx];}
					else{next=-1;}
					//				System.out.println("<A2>\ti="+i+", lidx="+lidx+", lims.length="+lims.length+", next="+next+", sarray.length="+sarray.length+", slen="+slen);
					//				System.out.prontln(sb2);
				}
				//			System.out.println("<A3>\ti="+i+", lidx="+lidx+", lims.length="+lims.length+", next="+next+", sarray.length="+sarray.length+", slen="+slen);

				if(i==numOverCutoffPlusOne){
					numOver50=ssum;
					basesOver50=slen;
					fractionOver50=slen*100f/scaflen;
				}

				//			long a=carray[i];
				long b=sarray[i];
				long c=scarray1[i];
				long d=scarray2[i];

				if(b>0){
					csum-=c;
					ssum-=b;
					clen-=d;
					slen-=(b*i);
				}

				if(ln50==-1 && slen<=s50){
					ln50=i;
					ll50=ssum+b;
				}

				if(ln90==-1 && slen<=s90){
					ln90=i;
					ll90=ssum+b;
				}
				//			System.out.println("<A4>\tb="+b+", c="+c+", d="+d+", csum="+csum+", ssum="+ssum+", clen="+clen+", slen="+slen);

				//			System.out.println("<A5>\ti="+i+", lidx="+lidx+", lims.length="+lims.length+", next="+next+", sarray.length="+sarray.length+", slen="+slen);

			}

			for(Triple tp : tlist){
				//			assert(false) : tlist;
				while(tp.length>=next && lidx<lims.length){
					
					prevLine=formatX(next, ssum, csum, slen, clen);
//					if(!skipDuplicateLines || ssum!=prevSsum || csum!=prevCsum || slen!=prevSlen || clen!=prevClen){
						if(prevLine!=null && (next==0 || next>=minScaffold)){
							sb2.append(prevLine);
							sb2.append('\n');
							prevLine=null;
						}
//					}
					prevSsum=ssum; prevCsum=csum; prevSlen=slen; prevClen=clen;
					
					lidx++;
					if(lidx<lims.length){next=lims[lidx];}
					//				System.out.println("<B>\n"+sb2+"\ni="+"?"+", lidx="+lidx+", lims.length="+lims.length+", next="+next);
					//				else{next=-1;}
				}

				if(numOver50==0 && tp.length>numOverCutoff){
					numOver50=ssum;
					basesOver50=slen;
					fractionOver50=slen*100f/scaflen;
				}

				//			long a=carray[i];
				long b=tp.length;
				long c=tp.contigs;
				long d=tp.contiglen;

				if(b>0){
					csum-=c;
					ssum-=1;
					clen-=d;
					slen-=b;
				}

				if(ln50==-1 && slen<=s50){
					ln50=b;
					ll50=ssum+1;
				}

				if(ln90==-1 && slen<=s90){
					ln90=b;
					ll90=ssum+1;
				}

			}
			if(prevLine!=null){
				sb2.append(prevLine);
				prevLine=null;
			}
		}
		
		clen=contiglen;
		csum=contigs;
		for(int i=0; i<carray.length && clen>0; i++){
			long a=carray[i];
			
			csum-=a;
			clen-=a*i;
			
			if(cn50==-1 && clen<=c50){
				cn50=i;
				cl50=csum+a;
			}
			
			if(cn90==-1 && clen<=c90){
				cn90=i;
				cl90=csum+a;
			}
		}

		for(int i=0; i<llist.size && clen>0; i++){
			long a=larray[i];
			
			csum-=1;
			clen-=a;
			
			if(cn50==-1 && clen<=c50){
				cn50=a;
				cl50=csum+1;
			}
			
			if(cn90==-1 && clen<=c90){
				cn90=a;
				cl90=csum+1;
			}
		}

		ln50=Tools.max(ln50, 0);
		ll50=Tools.max(ll50, 0);
		cn50=Tools.max(cn50, 0);
		cl50=Tools.max(cl50, 0);
		
		ln90=Tools.max(ln90, 0);
		ll90=Tools.max(ll90, 0);
		cn90=Tools.max(cn90, 0);
		cl90=Tools.max(cl90, 0);

		double sLogSum=calcLogSumScaffolds(slist, tlist, logSumOffset, 2);
		double cLogSum=calcLogSumContigs(clist, llist, logSumOffset, 2);
		double sPowSum=calcPowerSumScaffolds(slist, tlist, logSumOffset, powSumPower);
		double cPowSum=calcPowerSumContigs(clist, llist, logSumOffset, powSumPower);
		double sAssemblyScore=calcAssemblyScore(slist, tlist, assemblyScoreMinLen, assemblyScoreMaxLen, alignedFraction);
		
//		ByteStreamWriter tsw=new ByteStreamWriter((out==null ? "stdout" : out) , overwrite, append, false);
		TextStreamWriter tsw=new TextStreamWriter((out==null ? "stdout" : out) , overwrite, append, false);
		tsw.start();
		
		lastL50=ln50;
		lastL90=ln90;
		lastSize=contiglen;
		lastContigs=contigs;
		lastMaxContig=maxContig;
		lastLogSum=sLogSum;

		String sLogSumString=String.format((sLogSum/logSumOffset>=10000 ? "%.0f" : "%.3f"), sLogSum/logSumOffset);
		String cLogSumString=String.format((cLogSum/logSumOffset>=10000 ? "%.0f" : "%.3f"), cLogSum/logSumOffset);
		String sPowSumString=String.format((sPowSum/logSumOffset>=10000 ? "%.0f" : "%.3f"), sPowSum/logSumOffset);
		String cPowSumString=String.format((cPowSum/logSumOffset>=10000 ? "%.0f" : "%.3f"), cPowSum/logSumOffset);
		String sAssemblyScoreString=String.format((sAssemblyScore>=10000 ? "%.0f" : "%.3f"), sAssemblyScore);
		
		if(FORMAT<1){
			//Do nothing
		}else if(FORMAT<3){
			
			if(addfilename){sb.append("Filename:                           \t"+name+"\n");}
			sb.append("Main genome scaffold N/L50:         \t"+ll50+"/"+formatKB(ln50, 3, 0)+"\n");
			sb.append("Main genome contig N/L50:           \t"+cl50+"/"+formatKB(cn50, 3, 0)+"\n");
			if(printL90){
				sb.append("Main genome scaffold N/L90:         \t"+ll90+"/"+formatKB(ln90, 3, 0)+"\n");
				sb.append("Main genome contig N/L90:           \t"+cl90+"/"+formatKB(cn90, 3, 0)+"\n");
			}
			if(printExtended){
				sb.append("Main genome scaffold logsum:        \t"+sLogSumString+"\n");
				sb.append("Main genome scaffold powersum:      \t"+sPowSumString+"\n");
			}
			if(printExtended || printAssemblyScore){
				sb.append("Main genome assembly score:         \t"+sAssemblyScoreString+"\n");
			}
			sb.append("Max scaffold length:                \t"+formatKB(maxScaf, 3, 0)+"\n");
			sb.append("Max contig length:                  \t"+formatKB(maxContig, 3, 0)+"\n");
			sb.append("Number of scaffolds > 50 KB:        \t"+numOver50+"\n");
			if(printExtended){
				sb.append("Number of bases in scaffolds>50 KB: \t"+basesOver50+"\n");
			}
			
			sb.append("% main genome in scaffolds > 50 KB: \t"+Tools.format("%.2f%%", fractionOver50)+"\n");
			if(printheadersize){sb.append("Header:\t"+formatKB(HEADERLENSUM, 3, 0)+(HEADERLENSUM<1000 ? " bytes" : ""));}
			
			//		System.out.println();
			//		System.out.println("Scaffolds: "+shared.Vector.sum(slist.array));
			//		for(int i=0; i<slist.size; i++){
			//			if(slist.array[i]>0){System.out.print(i+":"+slist.array[i]+", ");}
			//		}
			//		System.out.println();
			//		System.out.println("Contigs:"+shared.Vector.sum(clist.array));
			//		for(int i=0; i<clist.size; i++){
			//			if(clist.array[i]>0){System.out.print(i+":"+clist.array[i]+", ");}
			//		}
			
			if(GCFORMAT==0){
				//Print nothing
			}else{
				if(GCFORMAT==1 || GCFORMAT==3 || GCFORMAT==4 || GCFORMAT==5){
					tsw.println("A\tC\tG\tT\tN\tIUPAC\tOther\tGC\tGC_stdev");
				}else{
					tsw.println("GC\tGC_stdev");
				}
				tsw.println(toString3(new StringBuilder(/*"Base Content"*/), counts, gc_std));
			}
			
			tsw.println(sb);
			tsw.println(sb2);
		}else if(FORMAT==3 || FORMAT==6){
			
			if(useheader){
				if(FORMAT==6){sb.append('#');}
				if(N_UNDERSCORE){sb.append("n_");}
				sb.append("scaffolds\t");
				if(N_UNDERSCORE){sb.append("n_");}
				sb.append("contigs\t");
				sb.append("scaf_bp\t");
				sb.append("contig_bp\t");
				sb.append("gap_pct\t");
				sb.append("scaf_N50\t");
				sb.append("scaf_L50\t");
				sb.append("ctg_N50\t");
				sb.append("ctg_L50\t");
				if(printL90){
					sb.append("scaf_N90\t");
					sb.append("scaf_L90\t");
					sb.append("ctg_N90\t");
					sb.append("ctg_L90\t");
				}
				if(printExtended){
					sb.append("scaf_logsum\t");
					sb.append("scaf_powsum\t");
					sb.append("ctg_logsum\t");
					sb.append("ctg_powsum\t");
				}
				if(printExtended || printAssemblyScore){
					sb.append("asm_score\t");
				}
				sb.append("scaf_max\t");
				sb.append("ctg_max\t");
				sb.append("scaf_n_gt50K\t");
				if(printExtended){
					sb.append("scaf_l_gt50k\t");
				}
				sb.append("scaf_pct_gt50K\t");
				sb.append("gc_avg\t");
				sb.append("gc_std");
				if(addfilename){sb.append("\tfilename");}
				
				sb.append("\n");
			}
			
			sb.append(scaffolds+"\t");
			sb.append(contigs+"\t");
			sb.append(scaflen+"\t");
			sb.append(Tools.format("%d",contiglen)+"\t");
			sb.append(Tools.format("%.3f",(scaflen-contiglen)*100f/scaflen)+"\t");
			sb.append(ll50+"\t");
			sb.append(formatKB(ln50, 3, 0)+"\t");
			sb.append(cl50+"\t");
			sb.append(formatKB(cn50, 3, 0)+"\t");
			if(printL90){
				sb.append(ll90+"\t");
				sb.append(formatKB(ln90, 3, 0)+"\t");
				sb.append(cl90+"\t");
				sb.append(formatKB(cn90, 3, 0)+"\t");
			}
			if(printExtended){
				sb.append(sLogSumString+"\t");
				sb.append(sPowSumString+"\t");
				sb.append(cLogSumString+"\t");
				sb.append(cPowSumString+"\t");
			}
			if(printExtended || printAssemblyScore){
				sb.append(sAssemblyScoreString+"\t");
			}
			sb.append(formatKB(maxScaf, 3, 0)+"\t");
			sb.append(formatKB(maxContig, 3, 0)+"\t");
			sb.append(numOver50+"\t");
			if(printExtended){
				sb.append(basesOver50+"\t");
			}
			sb.append(Tools.format("%.3f", fractionOver50)+"\t");
			sb.append(Tools.format("%.5f", (counts[1]+counts[2])*1.0/(counts[0]+counts[1]+counts[2]+counts[3]))+"\t");
			sb.append(Tools.format("%.5f", gc_std));
			if(addfilename){sb.append('\t').append(name);}
			
			tsw.println(sb);
		}else if(FORMAT==8){//JSON
			
			JsonObject jo=new JsonObject();
			jo.add("scaffolds", scaffolds);
			jo.add("contigs", contigs);
			jo.add("scaf_bp", scaflen);
			jo.add("contig_bp", contiglen);
			jo.addLiteral("gap_pct", (scaflen-contiglen)*100f/scaflen, 5);
			jo.add("scaf_N50", ll50);
			jo.add("scaf_L50", ln50);
			jo.add("ctg_N50", cl50);
			jo.add("ctg_L50", cn50);

			jo.add("scaf_N90", ll90);
			jo.add("scaf_L90", ln90);
			jo.add("ctg_N90", cl90);
			jo.add("ctg_L90", cn90);

			jo.addLiteral("scaf_logsum", sLogSumString);
			jo.addLiteral("scaf_powsum", sPowSumString);
			jo.addLiteral("ctg_logsum", cLogSumString);
			jo.addLiteral("ctg_powsum", cPowSumString);
			jo.addLiteral("asm_score", sAssemblyScoreString);
			jo.add("scaf_max", maxScaf);
			jo.add("ctg_max", maxContig);
			jo.add("scaf_n_gt50K", numOver50);
			jo.add("scaf_l_gt50k", basesOver50);

			jo.add("scaf_pct_gt50K", fractionOver50);
			jo.addLiteral("gc_avg", (counts[1]+counts[2])*1.0/(counts[0]+counts[1]+counts[2]+counts[3]), 5);
			jo.addLiteral("gc_std", gc_std, 5);
			jo.add("filename", name);
			
			tsw.println(jo.toString());
		}else if(FORMAT==4){
			
			if(useheader){

			if(N_UNDERSCORE){sb.append("n_");}
			sb.append("scaffolds\t");
//			sb.append("contigs\t");
			sb.append("scaf_bp\t");
//			sb.append("contig_bp\t");
//			sb.append("gap_pct\t");
			sb.append("scaf_N50\t");
			sb.append("scaf_L50\t");
			if(printL90){
				sb.append("scaf_N90\t");
				sb.append("scaf_L90\t");
			}
			if(printExtended){
				sb.append("scaf_logsum\t");
				sb.append("scaf_powsum\t");
			}

			if(printExtended || printAssemblyScore){
				sb.append("asm_score\t");
			}
//			sb.append("ctg_N50\t");
//			sb.append("ctg_L50\t");
			sb.append("scaf_max\t");
//			sb.append("ctg_max\t");
			sb.append("scaf_n_gt50K");
			if(printExtended){
				sb.append("\tscaf_l_gt50k");
			}
			sb.append("\tscaf_pct_gt50K");
			sb.append("\tgc_avg");
			sb.append("\tgc_std");
			if(addfilename){sb.append("\tfilename");}

			sb.append("\n");
			}
			
			sb.append(scaffolds+"\t");
//			sb.append(contigs+"\t");
			sb.append(scaflen+"\t");
//			sb.append(Tools.format("%d",contiglen)+"\t");
//			sb.append(Tools.format("%.3f",(scaflen-contiglen)*100f/scaflen)+"\t");
			sb.append(ll50+"\t");
			sb.append(formatKB(ln50, 3, 0)+"\t");
			if(printL90){
				sb.append(ll90+"\t");
				sb.append(formatKB(ln90, 3, 0)+"\t");
			}
			if(printExtended){
				sb.append(sLogSumString+"\t");
				sb.append(sPowSumString+"\t");
			}

			if(printExtended || printAssemblyScore){
				sb.append(sAssemblyScoreString+"\t");
			}
//			sb.append(cl50+"\t");
//			sb.append(formatKB(cn50, 3, 0)+"\t");
			sb.append(formatKB(maxScaf, 3, 0)+"\t");
//			sb.append(formatKB(maxContig, 3, 0)+"\t");
			sb.append(numOver50);
			if(printExtended){
				sb.append('\t').append(basesOver50);
			}
			sb.append(Tools.format("\t%.3f", fractionOver50));

			sb.append(Tools.format("\t%.5f", (counts[1]+counts[2])*1.0/(counts[0]+counts[1]+counts[2]+counts[3]))+"\t");
			sb.append(Tools.format("\t%.5f", gc_std));
			
			if(addfilename){sb.append('\t').append(name);}
			tsw.println(sb);
		}else if(FORMAT==5){
			
			if(useheader){
//			sb.append("scaffolds\t");
			if(N_UNDERSCORE){sb.append("n_");}
			sb.append("contigs\t");
//			sb.append("scaf_bp\t");
			sb.append("contig_bp\t");
			sb.append("gap_pct\t");
//			sb.append("scaf_N50\t");
//			sb.append("scaf_L50\t");
			sb.append("ctg_N50\t");
			sb.append("ctg_L50\t");
			if(printL90){
				sb.append("ctg_N90\t");
				sb.append("ctg_L90\t");
			}
			if(printExtended){				
				sb.append("ctg_logsum\t");
				sb.append("ctg_powsum\t");
			}

			if(printExtended || printAssemblyScore){
				sb.append("asm_score\t");
			}
//			sb.append("scaf_max\t");
			sb.append("ctg_max\t");
//			sb.append("scaf_n_gt50K\t");
//			sb.append("scaf_pct_gt50K\t");
			sb.append("gc_avg\t");
			sb.append("gc_std");
			if(addfilename){sb.append("\tfilename");}

			sb.append("\n");
			}
			
//			sb.append(scaffolds+"\t");
			sb.append(contigs+"\t");
//			sb.append(scaflen+"\t");
			sb.append(Tools.format("%d",contiglen)+"\t");
			sb.append(Tools.format("%.3f",(scaflen-contiglen)*100f/scaflen)+"\t");
//			sb.append(ll50+"\t");
//			sb.append(formatKB(ln50, 3, 0)+"\t");
			sb.append(cl50+"\t");
			sb.append(formatKB(cn50, 3, 0)+"\t");
			if(printL90){
				sb.append(cl90+"\t");
				sb.append(formatKB(cn90, 3, 0)+"\t");
			}
			if(printExtended){
				sb.append(cLogSumString+"\t");
				sb.append(cPowSumString+"\t");
			}
			if(printExtended || printAssemblyScore){
				sb.append(sAssemblyScoreString+"\t");
			}
//			sb.append(formatKB(maxScaf, 3, 0)+"\t");
			sb.append(formatKB(maxContig, 3, 0)+"\t");
//			sb.append(numOver50+"\t");
//			sb.append(Tools.format("%.3f", fractionOver50)+"\t");
			sb.append(Tools.format("%.5f", (counts[1]+counts[2])*1.0/(counts[0]+counts[1]+counts[2]+counts[3]))+"\t");
			sb.append(Tools.format("%.5f", gc_std));
			if(addfilename){sb.append('\t').append(name);}
			tsw.println(sb);
		}else if(FORMAT==7){

			if(addfilename){sb.append("Filename:                           \t"+name+"\n");}
			sb.append("Main genome contig N/L50:           \t"+cl50+"/"+formatKB(cn50, 3, 0)+"\n");
			if(printL90){
				sb.append("Main genome contig N/L90:           \t"+cl90+"/"+formatKB(cn90, 3, 0)+"\n");
			}
			if(printExtended){
//				sb.append("Main genome contig logsum:          \t"+cLogSumString+"\n");
			}
			sb.append("Max contig length:                  \t"+formatKB(maxContig, 3, 0)+"\n");
			sb.append("Number of contigs > 1 KB:           \t"+numOver50+"\n");
			sb.append("% main genome in contigs > 1 KB:    \t"+Tools.format("%.2f%%", fractionOver50)+"\n");
			if(printheadersize){sb.append("Header:\t"+formatKB(HEADERLENSUM, 3, 0)+(HEADERLENSUM<1000 ? " bytes" : ""));}
			
			if(GCFORMAT==0){
				//Print nothing
			}else{
				if(GCFORMAT==1 || GCFORMAT==3 || GCFORMAT==4 || GCFORMAT==5){
					tsw.println("A\tC\tG\tT\tGC\tGC_stdev");
				}else{
					tsw.println("GC\tGC_stdev");
				}
				tsw.println(toString3(new StringBuilder(/*"Base Content"*/), counts, gc_std));
			}

			tsw.println(sb);
			tsw.println(sb2);
		}
		tsw.poisonAndWait();
		
		if(showspeed){
			if(!printheadersize){outstream.println("Header:\t"+formatKB(HEADERLENSUM, 3, 0)+(HEADERLENSUM<1000 ? " bytes" : ""));}
			outstream.println("Time: \t"+t);
			long bytes=new File(in).length();
			if(bytes<1){bytes=LIMSUM;}
			double mbps=bytes*1000d/t.elapsed;
			double mbpps=shared.Vector.sum(counts)*1000d/t.elapsed;
			outstream.println(Tools.format("Speed:\t%.2f MBytes/s",mbps));
			outstream.println(Tools.format("      \t%.2f MBases/s",mbpps));
		}
		
		if(bbmapkmer>0){
			outstream.println("BBMap minimum memory estimate at k="+bbmapkmer+":     "+estimateBBMapMemory(counts, scaffolds, HEADERLENSUM, bbmapkmer));
		}
		
		MetadataWriter.write(null, scaffolds, scaflen, 0, 0, false);
	}
	
	private static long bbmapMemoryBytes(long[] acgtn, long scaffolds,
			long headerlen, int k) {

		long keyspace=(1L<<(2*k));
		long defined=acgtn[0]+acgtn[1]+acgtn[2]+acgtn[3];
		long undefined=acgtn[4];
		long midpad=(scaffolds*(FastaToChromArrays2.MID_PADDING));
		long total=defined+undefined+midpad;
		int chromlen=FastaToChromArrays2.MAX_LENGTH-FastaToChromArrays2.END_PADDING-FastaToChromArrays2.START_PADDING;
		int chroms=(int)(total/chromlen);
		int chromsperblock=Integer.MAX_VALUE/chromlen;
		int blocks=(chroms+chromsperblock-1)/chromsperblock;
		long memperblock=keyspace*4;
		long memforcounts=keyspace*4;
		
		long mem=0;
		mem+=total; //reference bulk, including inter-scaffold padding
		mem+=(chroms*(FastaToChromArrays2.END_PADDING+FastaToChromArrays2.START_PADDING)); //reference tip padding
		mem+=headerlen; //Header name byte arrays
		mem+=(scaffolds*(4+4+4+16+8)); //Other structures for scaffold info
		mem+=(blocks*(memperblock)); //start array for each block
		mem+=memforcounts; //count array
		mem+=(defined*4); //key lists
		mem=(long)(mem/0.66); //Expand to compensate for garbage collection
		if(k>13){mem=mem+1000000000;}
		return mem;
	}
	
	private static CharSequence estimateBBMapMemory(long[] acgtn, long scaffolds,
			long headerlen, int k) {
		long mem=180+bbmapMemoryBytes(acgtn, scaffolds, headerlen, k)/1000000; //in megabytes
		if(mem>4000){
			return "-Xmx"+((mem+1500)/1000)+"g \t(at least "+(long)Math.ceil((((mem+1500)/0.85)/1000))+" GB physical RAM)";
		}else if(mem>2100){
			return "-Xmx"+((mem+999)/1000)+"g \t(at least "+(long)Math.ceil((((mem+999)/0.85)/1000))+" GB physical RAM)";
		}else{
			return "-Xmx"+(((((mem*11)/8+50))/10)*10)+"m \t(at least "+((((long)(((mem*10)/8+50)/0.82))/10)*10)+" MB physical RAM)";
		}
	}


	public static long bench(InputStream is) throws IOException{
		final byte[] buf=new byte[32768];
		long sum=0;
		for(long len=is.read(buf); len>0; len=is.read(buf)){sum+=len;}
		return sum;
	}
	
	
	
	private static void writeHistFile(String fname, LongList slist, ArrayList<Triple> tlist, boolean ascending){
		if(fname==null){return;}
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		tsw.print("#scaffolds\tlength\n");
		long num=0, len=0;
		
		final long[] array=slist.array;
		final int lim=slist.size;
		
		StringBuilder sb=new StringBuilder(32);
		
		if(ascending){
			for(int i=0; i<lim; i++){
				final long a=array[i];
				if(a>0){
					len+=(a*i);
					num+=a;
					sb.append(num);
					sb.append('\t');
					sb.append(len);
					sb.append('\n');
					tsw.print(sb.toString());
					sb.setLength(0);
				}
			}

			if(tlist!=null){
				for(Triple t : tlist){
					len+=t.length;
					num++;
					sb.append(num);
					sb.append('\t');
					sb.append(len);
					sb.append('\n');
					tsw.print(sb.toString());
					sb.setLength(0);
				}
			}
		}else{

			if(tlist!=null){
				for(int i=tlist.size()-1; i>=0; i--){
					Triple t=tlist.get(i);
					len+=t.length;
					num++;
					sb.append(num);
					sb.append('\t');
					sb.append(len);
					sb.append('\n');
					tsw.print(sb.toString());
					sb.setLength(0);
				}
			}
			for(int i=lim-1; i>=0; i--){
				final long a=array[i];
				if(a>0){
					len+=(a*i);
					num+=a;
					sb.append(num);
					sb.append('\t');
					sb.append(len);
					sb.append('\n');
					tsw.print(sb.toString());
					sb.setLength(0);
				}
			}
		}
		tsw.poisonAndWait();
	}
	
	private static String toString2(StringBuilder sb, long[] counts){
		final long a=counts[0], c=counts[1], g=counts[2], t=counts[3], iupac=counts[4], n=counts[5], other=counts[6], control=counts[7];
		long sumDef=a+c+g+t;
		long sumAll=sumDef+iupac+n+other;
		double invDef=1.0/sumDef, invAll=1.0/sumAll;
		double iupacD=iupac*invAll;
		double otherD=other*invAll;
		if(iupac>0 && iupacD<0.0001){iupacD=0.0001;}
		if(other>0 && otherD<0.0001){otherD=0.0001;}
		if(GCFORMAT==0 || GCFORMAT==1){
			return sb.append(Tools.format("\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
					sumAll, a*invDef, c*invDef, g*invDef, t*invDef, n*invAll, iupacD, otherD, (g+c)*invDef)).toString();
		}else if(GCFORMAT==2){
			return sb.append(Tools.format("\t%.4f\n", (g+c)*invDef)).toString();
		}else if(GCFORMAT==3){
			return sb.append(Tools.format("\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
					sumAll, a*invDef, c*invDef, g*invDef, t*invDef, n*invAll, iupacD, otherD, (g+c)*invDef)).toString();
		}else if(GCFORMAT==4){
			return sb.append(Tools.format("\t%d\t%.4f\n", sumAll, (g+c)*invDef)).toString();
		}else if(GCFORMAT==5){
			return sb.append(Tools.format("\t%d\t%.4f\t%.2f\t%.2f\n", sumAll, (g+c)*invDef, logsum(sumAll), powersum(sumAll))).toString();
		}else{
			throw new RuntimeException("Unknown format.");
		}
	}
	
	private static String toString3(StringBuilder sb, long[] counts, double gc_std){
		final long a=counts[0], c=counts[1], g=counts[2], t=counts[3], iupac=counts[4], n=counts[5], other=counts[6], control=counts[7];
		long sumDef=a+c+g+t;
		long sumAll=sumDef+iupac+n+other;
		double invDef=1.0/sumDef, invAll=1.0/sumAll;
		double iupacD=iupac*invAll;
		double otherD=other*invAll;
		if(iupac>0 && iupacD<0.0001){iupacD=0.0001;}
		if(other>0 && otherD<0.0001){otherD=0.0001;}
		if(FORMAT==7){
			return sb.append(Tools.format("%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
					a*invDef, c*invDef, g*invDef, t*invDef, (g+c)*invDef, gc_std)).toString();
		}else if(GCFORMAT==0 || GCFORMAT==1 || GCFORMAT==3 || GCFORMAT==4 || GCFORMAT==5){
			if(GC_FRACTION){
			return sb.append(Tools.format("%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
					a*invDef, c*invDef, g*invDef, t*invDef, n*invAll, iupacD, otherD, (g+c)*invDef, gc_std)).toString();
			}
			return sb.append(Tools.format("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\n",
					a, c, g, t, n, iupac, other, (g+c)*invDef, gc_std)).toString();
		}else if(GCFORMAT==2){
			return sb.append(Tools.format("%.4f\t%.4f\n", (g+c)*invDef, gc_std)).toString();
		}else{
			throw new RuntimeException("Unknown format.");
		}
		
	}
	
	
	private static final String formatX(int next, long ssum, long csum, long slen, long clen){
		float cov=clen*100f/slen;
		
		final String s;
		if(FORMAT<=1){
			s=formatKB_all(next, 1, 7)+" \t"+formatComma(ssum, 14)+"\t"+formatComma(csum, 14)+"\t"+formatComma(slen, 14)+"\t"+formatComma(clen, 14)+"\t"+formatPercent(cov);
		}else if(FORMAT==7){
			if(next==0){return null;}
			s=formatKB_all(next, 1, 7)+" \t"+formatComma(csum, 14)+"\t"+formatComma(clen, 14);
		}else if(FORMAT>=2){
			s=formatKB_all(next, 1, 7)+" \t"+formatComma(ssum, 14)+"\t"+formatComma(csum, 14)+"\t"+formatComma(slen, 14)+"\t"+formatComma(clen, 14)+"\t"+formatPercent(cov);
		}else{
			throw new RuntimeException("Unknown format: "+FORMAT);
		}
		return s;
	}
	
	private static final String formatKB(long x, int precision, int width){
		String s;
		if(FORMAT>=2 || x<1000){
			s=Long.toString(x);
		}else if(x<1000000){
			s=Tools.format("%."+precision+"f",x/1000f);
			while(s.contains(".") && (s.endsWith("0") || s.endsWith("."))){s=s.substring(0, s.length()-1);}
			s=s+" Kbp";
		}else{
			s=Tools.format("%."+precision+"f",x/1000000f);
			while(s.contains(".") && (s.endsWith("0") || s.endsWith("."))){s=s.substring(0, s.length()-1);}
			s=s+" Mbp";
		}
		while(s.length()<width){s=" "+s;}
		return s;
	}
	
	private static final String formatKB_all(long x, int precision, int width){
		String s;
			
		if(x==0){s="All";}
		else if(FORMAT>=2 || x<1000){
			s=Long.toString(x);
		}else if(x<1000000){
			s=Tools.format("%."+precision+"f",x/1000f);
			while(s.contains(".") && (s.endsWith("0") || s.endsWith("."))){s=s.substring(0, s.length()-1);}
			s=s+" Kbp";
		}else{
			s=Tools.format("%."+precision+"f",x/1000000f);
			while(s.contains(".") && (s.endsWith("0") || s.endsWith("."))){s=s.substring(0, s.length()-1);}
			s=s+" Mbp";
		}
		
		while(s.length()<width){s=" "+s;}
		return s;
	}
	
	private static final StringBuilder formatComma(long x, int width){
		StringBuilder sb=new StringBuilder(width);
		if(FORMAT<=1){
			sb.append(x%1000);
			x/=1000;
			int len=3;
			while(x>0){
				while(sb.length()<len){sb.insert(0, '0');}
				sb.insert(0, (x%1000)+",");
				x/=1000;
				len+=4;
			}
		}else if(FORMAT>=2){
			sb.append(x);
		}else{
			throw new RuntimeException("Unknown format: "+FORMAT);
		}
		while(sb.length()<width){
			sb.insert(0, ' ');
		}
		return sb;
	}
	
	private static final String formatPercent(float x){
		String s=Tools.format("%.2f%%", x);
		while(s.length()<8){s=" "+s;}
		return s;
	}
	
	protected void reset(){
//		clist=null;
//		slist=null;
//		sclist1=null;
//		sclist2=null;
//
//		gchistArray=null;
//		gc_std=0;
//
//		llist=null;
//		tlist=null;

		clist.clear();
		slist.clear();
		sclist1.clear();
		sclist2.clear();
		llist.clear();
		tlist.clear();

		Arrays.fill(gchistArray, 0);
		Arrays.fill(gchist_by_base, 0);
		
		gchistArray_downsampled=null;
		gchist_by_base_downsampled=null;
		
		gc_std=0;
		gc_bb_std=0;
		
		LIMSUM=0;
		HEADERLENSUM=0;
	}
	
	/**
	 * @return charToNum array
	 */
	public static final byte[] makeCharToNum() {
		byte[] r=new byte[256];
		Arrays.fill(r, (byte)6);
		r['a']=r['A']=0;
		r['c']=r['C']=1;
		r['g']=r['G']=2;
		r['t']=r['T']=3;
		r['u']=r['U']=3;
		r['n']=r['N']=5;
		r['x']=r['X']=4;
		for(byte b : AminoAcid.degenerateBases){
			if(b!=' '){
				r[b]=r[Tools.toLowerCase(b)]=4;
			}
		}
		r['\n']=r['\r']=r['>']=r['@']=r['+']=7;
		return r;
	}
	
	/*--------------------------------------------------------------*/
	
	private static final byte[] charToNum=makeCharToNum();
	public static int GCFORMAT=1;
	public static boolean GC_FRACTION=true;
	public static int FORMAT=1;
	private static long cutoff=1000000;
	
	private static long LIMSUM=0;
	private static long HEADERLENSUM=0;
	
	private static int bbmapkmer=0;//13;
	public static boolean overwrite=true;
	public static boolean append=false;
	public static boolean useheader=true;
	public static boolean addfilename=false;
	public static boolean printL90=true;
	public static boolean printExtended=false;
	public static int logSumOffset=1000;
	public static double powSumPower=0.25;
	public static double logSumBase=2;
	public static boolean squareLog=false;
	public static double logPower=1.0;
	public static boolean showspeed=false;//true;
	public static boolean printheadersize=false;
	public static boolean skipDuplicateLines=true;
	public static boolean N_UNDERSCORE=true;

	public static long assemblyScoreMinLen=2000;
	public static long assemblyScoreMaxLen=50000;
	public static double alignedFraction=0;
	public static boolean printAssemblyScore=false;

	private static final byte slashr='\r', slashn='\n', carrot='>', at='@', noref='N', noref2='n';
	
	/*--------------------------------------------------------------*/
	
	private boolean benchmark=false;
	private String in=null, out=null, gc=null, gchistFile=null, scaffoldHistFile=null;
	private int maxNs=-1;
	
	/** Number of decimal places for GC histogram */
	private final int gchistDecimals1;
	
	/** Number of bins for output (subsampled) GC content histogram */
	private final int gcbins;
	
	/** Number of bins for internal GC content histogram */
	private final int gcbins2;
	
	/** Minimum scaffold length to count */
	private final int minScaffold;
	
	/** Number of contigs of length x */
	private final LongList clist;
	
	/** Number of scaffolds of length x */
	private final LongList slist;
	
	/** Sum of contigs per scaffold of length x */
	private final LongList sclist1;
	
	/** Sum of contig lengths per scaffold of length x */
	private final LongList sclist2;
	
	/** List of contig lengths for contigs at least cutoff in length */
	private final LongList llist;
	
	/** List of scaf len, contigs, contig sum for scaffolds at least cutoff in length */
	private final ArrayList<Triple> tlist;

	/** Downsampled gc histogram */
	private final long[] gchistArray;

	/** Downsampled gc histogram */
	private long[] gchistArray_downsampled;
	
	/** gc standard deviation */
	private double gc_std;
	
	/** Downsampled gc histogram, using base counts rather than scaffold counts */
	private final long[] gchist_by_base;
	
	/** Downsampled gc histogram, using base counts rather than scaffold counts */
	private long[] gchist_by_base_downsampled;
	
	/** gc standard deviation, using base counts rather than scaffold counts */
	private double gc_bb_std;

	public static double lastLogSum;
	public static long lastL90;
	public static long lastL50;
	public static long lastSize;
	public static long lastContigs;
	public static long lastMaxContig;
	
	public static PrintStream outstream=System.out;
	
	/*--------------------------------------------------------------*/
	
	private static class Triple implements Comparable<Triple>{
		
		Triple(long len_, long contigs_, long contiglen_){
			length=len_;
			contigs=contigs_;
			contiglen=contiglen_;
		}
		
		@Override
		public int compareTo(Triple o) {
			if(length>o.length){return 1;}
			if(length<o.length){return -1;}
			return (int)(contiglen-o.contiglen);
		}
		
		@Override
		public boolean equals(Object o){return equals((Triple)o);}
		
		public boolean equals(Triple o){
			return length==o.length && contiglen==o.contiglen;
		}
		
		@Override
		public int hashCode() {
			assert(false) : "This class should not be hashed.";
			return super.hashCode();
		}
		
		@Override
		public String toString(){return length+","+contigs+","+contiglen;}

		public final long length;
		public final long contigs;
		public final long contiglen;
		
	}
	
	/*--------------------------------------------------------------*/
	

	
	/*
	 
	fasta_stats2.linux -n <number of N between contigs> contigs.fa
	e.g.

	fasta_stats2.linux -n 0 contigs.fa # for aplg assemblies
	fasta_stats2.linux -n 10 contigs.fa # for velvet
	  
	  
	 
	 Main genome scaffold total: 1610
	 Main genome contig total:   7844
	 Main genome scaffold sequence total: 726.6 MB
	 Main genome contig sequence total:   689.4 MB (->  5.1% gap)
	 Main genome scaffold N/L50: 6/62.2 MB
	 Main genome contig N/L50:   331/429.0 KB
	 Number of scaffolds > 50 KB: 122
	 % main genome in scaffolds > 50 KB: 98.9%

	  Minimum    Number    Number     Total        Total     Scaffold
	 Scaffold      of        of      Scaffold      Contig     Contig
	  Length   Scaffolds  Contigs     Length       Length    Coverage
	 --------  ---------  -------  -----------  -----------  --------
	     All     1,610      7,844  726,616,606  689,442,341    94.88%
	    1 kb     1,610      7,844  726,616,606  689,442,341    94.88%
	  2.5 kb     1,468      7,677  726,334,758  689,171,164    94.88%
	    5 kb       537      6,496  723,058,922  685,949,825    94.87%
	   10 kb       321      6,176  721,557,480  684,511,419    94.87%
	   25 kb       138      5,900  718,873,396  681,879,275    94.85%
	   50 kb       122      5,854  718,322,923  681,420,273    94.86%
	  100 kb        83      5,660  715,543,850  679,452,337    94.96%
	  250 kb        47      5,326  709,779,897  675,162,461    95.12%
	  500 kb        32      5,073  704,645,704  671,472,605    95.29%
	    1 mb        19      4,735  695,996,631  664,862,860    95.53%
	  2.5 mb        15      4,587  689,883,367  659,102,480    95.54%
	    5 mb        13      4,463  681,669,379  651,024,951    95.50%
	*/
	
}
