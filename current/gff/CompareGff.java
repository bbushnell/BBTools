package gff;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.ProkObject;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.StringNum;

/**
 * Compares GFF files for gene-calling accuracy assessment.
 * Evaluates genomic annotation precision by comparing reference GFF files against
 * query GFF files, calculating true positives, false positives, and false negatives
 * for feature starts and stops. Supports CDS, rRNA, and tRNA feature types with
 * strand-aware comparison and signal-to-noise ratio calculation.
 *
 * @author Brian Bushnell
 * @date October 3, 2018
 */
public class CompareGff {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Application entry point for GFF file comparison.
	 * @param args Command line arguments including reference and query file paths */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		CompareGff x=new CompareGff(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs CompareGff with command line arguments.
	 * Parses arguments, initializes file formats, and validates file accessibility.
	 * @param args Command line arguments for configuration
	 */
	public CompareGff(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;
			
			in=parser.in1;
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffin=FileFormat.testInput(in, FileFormat.GFF, null, true, true);
		ffref=FileFormat.testInput(ref, FileFormat.GFF, null, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses command line arguments into configuration settings.
	 * Supports ref=file, lines=count, verbose=boolean parameters.
	 * @param args Array of command line arguments
	 * @return Configured Parser object with parsed settings
	 */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("ref")){
				ref=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
//				ByteFile1.verbose=verbose;
//				ByteFile2.verbose=verbose;
//				ReadWrite.verbose=verbose;
			}else if(a.equals("ncrna") || a.equals("overlap")){
				ncrna=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(i==0 && arg.indexOf('=')<0){
				parser.in1=arg;
			}else if(i==1 && arg.indexOf('=')<0 && ref==null){
				ref=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		return parser;
	}
	
	/** Automatically adds or removes .gz/.bz2 extensions as needed.
	 * Validates that both input and reference files are specified. */
	private void fixExtensions(){
		in=Tools.fixExtension(in);
		ref=Tools.fixExtension(ref);
		if(in==null || ref==null){throw new RuntimeException("Error - at least two input files are required.");}
	}
	
	/** Validates that input and reference files exist and are readable.
	 * Throws RuntimeException if files cannot be accessed. */
	private void checkFileExistence(){
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(true, true, in, ref)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
	}
	
	/** Configures static ByteFile settings for optimal performance.
	 * Forces BF2 mode when using more than 2 threads for improved I/O. */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
//		if(!ByteFile.FORCE_MODE_BF2){
//			ByteFile.FORCE_MODE_BF2=false;
//			ByteFile.FORCE_MODE_BF1=true;
//		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin);
		
		processInner(bf);
		
		errorState|=bf.close();
		
		t.stop();
		
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Ref count:           \t"+refCount);
		outstream.println("Query count:         \t"+queryCount);

		outstream.println();
		outstream.println("Ref-relative counts:");
		outstream.println("True Positive Start: \t"+truePositiveStart+"\t"+(Tools.format("%.3f%%", truePositiveStart*100.0/refCount)));
		outstream.println("True Positive Stop:  \t"+truePositiveStop+"\t"+(Tools.format("%.3f%%", truePositiveStop*100.0/refCount)));
//		outstream.println("False Positive Start:\t"+falsePositiveStart+"\t"+(Tools.format("%.3f%%", falsePositiveStart*100.0/refCount)));
//		outstream.println("False Positive Stop: \t"+falsePositiveStop+"\t"+(Tools.format("%.3f%%", falsePositiveStop*100.0/refCount)));
		outstream.println("False Negative Start:\t"+falseNegativeStart+"\t"+(Tools.format("%.3f%%", falseNegativeStart*100.0/refCount)));
		outstream.println("False Negative Stop: \t"+falseNegativeStop+"\t"+(Tools.format("%.3f%%", falseNegativeStop*100.0/refCount)));

		outstream.println();
		outstream.println("Query-relative counts:");
		outstream.println("True Positive Start: \t"+truePositiveStart2+"\t"+(Tools.format("%.3f%%", truePositiveStart2*100.0/queryCount)));
		outstream.println("True Positive Stop:  \t"+truePositiveStop2+"\t"+(Tools.format("%.3f%%", truePositiveStop2*100.0/queryCount)));
		outstream.println("False Positive Start:\t"+falsePositiveStart2+"\t"+(Tools.format("%.3f%%", falsePositiveStart2*100.0/queryCount)));
		outstream.println("False Positive Stop: \t"+falsePositiveStop2+"\t"+(Tools.format("%.3f%%", falsePositiveStop2*100.0/queryCount)));
		
		outstream.println();
		outstream.println("SNR: \t"+Tools.format("%.4f", 10*Math.log10((truePositiveStart2+truePositiveStop2+0.1)/(falsePositiveStart2+falsePositiveStop2+0.1))));

		if(ncrna){processNcRNA();}

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@SuppressWarnings("unchecked")
	private void processInner(ByteFile bf){
		byte[] line=bf.nextLine();
		
		{
			ArrayList<GffLine> refLines=GffLine.loadGffFile(ffref, "CDS,rRNA,tRNA", true);

			refCount=refLines.size();
			lineMap=new HashMap<StringNum, GffLine>();
			startCountMap=new HashMap<StringNum, Integer>();
			stopCountMap=new HashMap<StringNum, Integer>();
			
			for(GffLine gline : refLines){
				final int stop=gline.trueStop();
				StringNum sn=new StringNum(gline.seqid, stop);
				lineMap.put(sn, gline);
				startCountMap.put(sn, 0);
				stopCountMap.put(sn, 0);
				assert(lineMap.get(sn)==gline);
//				assert(false) : "\n\nsn='"+sn+"'\n"+lineMap.containsKey(sn)+"\n"+lineMap.keySet();
			}
			if(verbose){
				System.err.println(lineMap);
				System.err.println(startCountMap);
				System.err.println(stopCountMap);
			}
		}

		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				
				final boolean valid=(line[0]!='#');
				if(valid){
					queryCount++;
					GffLine gline=new GffLine(line);
					processLine(gline);
				}
			}
			line=bf.nextLine();
		}
		
		for(Entry<StringNum, Integer> e : startCountMap.entrySet()){
			if(e.getValue()<1){
				falseNegativeStart++;
			}
		}
		for(Entry<StringNum, Integer> e : stopCountMap.entrySet()){
			if(e.getValue()<1){
				falseNegativeStop++;
			}
		}
	}
	
	private void processLine(GffLine gline){
//		boolean cds=gline.type.equals("CDS");
//		boolean trna=gline.type.equals("tRNA");
//		boolean rrna=gline.type.equals("rRNA");
//		if(!cds && !trna && !rrna){return;}
//		if(cds && !ProkObject.callCDS){return;}
//		if(trna && !ProkObject.calltRNA){return;}
//		if(rrna){
//			int type=gline.prokType();
//			if(ProkObject.processType(type)){return;}
//		}
		//prokType() is now null-safe (GffLine#001 fix): a truncated query line returns -1, processType(-1)=true -> not
		//skipped. Features are matched ref<->query by (seqid, trueStop): a stop match = TP-stop, +matching start = TP-start.
		int type=gline.prokType();
		if(!ProkObject.processType(type)){return;}

		final int stop=gline.trueStop();
		final int start=gline.trueStart();
		
//		System.err.println("Considering "+start+", "+stop);

		StringNum sn=new StringNum(gline.seqid, stop);
		GffLine refline=lineMap.get(sn);
		
		boolean fail=(refline==null || refline.strand!=gline.strand || !refline.type.equals(gline.type));
		if(fail){
			if(verbose){
				System.err.println("Can't find "+sn+"\n"+gline+"\n"+refline);
				assert(false) : "\n\nsn='"+sn+"'\n"+lineMap.containsKey(sn)+"\n"+lineMap.keySet();
			}
			falsePositiveStart++;
			falsePositiveStop++;
			falsePositiveStart2++;
			falsePositiveStop2++;
		}else{
			assert(stop==refline.trueStop());
			truePositiveStop++;
			truePositiveStop2++;
			stopCountMap.put(sn, stopCountMap.get(sn)+1);
			if(start==refline.trueStart()){
				truePositiveStart++;
				truePositiveStart2++;
				startCountMap.put(sn, startCountMap.get(sn)+1);
			}else{
				falsePositiveStart++;
				falsePositiveStart2++;
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      ncRNA Overlap Metrics   ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Overlap-based accuracy for ncRNA feature types (tRNA, rRNA), reported IN ADDITION to
	 * the exact-match counts above (which stay unchanged). A called feature is a true positive
	 * if it overlaps a reference feature of the same type, on the same seqid, by more than half
	 * of BOTH features' lengths (reciprocal 50% overlap == overlap &gt; 0.5*max(refLen, qryLen),
	 * using the eval_anticodon3 no-+1 length convention). For matched pairs it also reports
	 * 5'/3' boundary-exactness and offsets (strand-aware), and for tRNA the anticodon/amino-acid
	 * accuracy of the call by two paths -- model-name annotation and structural anticodon
	 * extraction -- against the reference product/Note annotation.
	 *
	 * Two deliberate deviations from the eval_anticodon3.sh reference:
	 * (1) RECIPROCAL overlap -- the script gates on the reference side only (overlap &gt; 0.5*refLen);
	 *     a one-sided gate lets a single genome-spanning prediction trivially match every reference,
	 *     so a general-purpose tool requires overlap on the query side too. The script gets away with
	 *     it because callgenes tRNA predictions are always ~70-120 bp, never genome-spanning.
	 * (2) SEQID-AWARE -- a call only matches a reference on the SAME contig, mandatory when the
	 *     query/ref hold many concatenated genomes; the per-genome python ignores seqid.
	 * Numbers match the reference for single-contig genomes whose calls are near reference length,
	 * and diverge only where a call is much longer/shorter than its reference (reciprocal rejects it)
	 * or on coincidental cross-contig overlap (seqid-awareness rejects it).
	 */
	private void processNcRNA(){
		overlapMetrics("tRNA", true);
		overlapMetrics("rRNA", false);
	}

	/** Computes and prints overlap-based metrics for one ncRNA feature type.
	 * @param type GFF feature type to score ("tRNA" or "rRNA")
	 * @param isTrna If true, also computes anticodon/amino-acid accuracy */
	private void overlapMetrics(String type, boolean isTrna){
		final ArrayList<GffLine> refList=GffLine.loadGffFile(ffref, type, false);
		final ArrayList<GffLine> qryList=GffLine.loadGffFile(ffin, type, false);

		//Index reference features per seqid, sorted by start, for fast overlap lookup
		final HashMap<String, ArrayList<GffLine>> refBySeqid=new HashMap<String, ArrayList<GffLine>>();
		for(GffLine r : refList){
			ArrayList<GffLine> list=refBySeqid.get(r.seqid);
			if(list==null){list=new ArrayList<GffLine>(); refBySeqid.put(r.seqid, list);}
			list.add(r);
		}
		for(ArrayList<GffLine> list : refBySeqid.values()){Collections.sort(list, START_COMP);}

		final int refCount=refList.size(), qryCount=qryList.size();
		long tp=0, fp=0;
		//Boundary tallies over matched pairs
		long bothExact=0, startExact=0, stopExact=0;
		long sumDs=0, sumDe=0, absDs=0, absDe=0;
		//Anticodon tallies (tRNA only)
		long extracted=0, skipUnk=0;
		long mAaC=0, mAaW=0, mAcC=0, mAcW=0;//model-name path: amino/anticodon correct/wrong
		long eAaC=0, eAaW=0, eAcC=0, eAcW=0;//structural-extraction path

		for(GffLine q : qryList){
			final int qs=q.start, qe=q.stop;
			String[] call=null;
			if(isTrna){
				call=calledAmino(q.attributes);
				if("1".equals(call[4])){extracted++;}//extraction coverage counts EVERY call, matched or not
			}

			//First reference on the same seqid overlapping by >50% of its own length wins
			GffLine match=null;
			final ArrayList<GffLine> candidates=refBySeqid.get(q.seqid);
			if(candidates!=null){
				for(GffLine r : candidates){
					if(r.start>qe){break;}//sorted by start: no later ref can overlap
					final int ov=Math.min(qe, r.stop)-Math.max(qs, r.start);//no-+1 length convention
					final int refLen=r.stop-r.start, qryLen=qe-qs;//no-+1 length convention
					//RECIPROCAL 50%: overlap must exceed half of BOTH features' lengths (== half of the
					//longer). This deviates from eval_anticodon3 (which gates on the reference side only)
					//per Brian: a one-sided gate lets one genome-spanning prediction trivially match every
					//reference. A general-purpose tool must require overlap on the query side too.
					if(2*ov>Math.max(refLen, qryLen)){match=r; break;}
				}
			}
			if(match==null){fp++; continue;}
			tp++;

			//Boundary accuracy (strand-aware: for minus strand, swap which end is 5'/3')
			int ds=qs-match.start, de=qe-match.stop;
			if(q.strand==GffLine.MINUS){final int tmp=ds; ds=-de; de=-tmp;}
			if(ds==0){startExact++;}
			if(de==0){stopExact++;}
			if(ds==0 && de==0){bothExact++;}
			sumDs+=ds; sumDe+=de; absDs+=Math.abs(ds); absDe+=Math.abs(de);

			if(isTrna){
				final String[] refAA=refAmino(match.attributes);//[aa, ac-or-null]
				final String raa=refAA[0], rac=refAA[1];
				if("UNK".equals(raa)){
					skipUnk++;
				}else{
					final String maa=call[0], mac=call[1], eaa=call[2], eac=call[3];
					if(aaMatch(maa, raa)){mAaC++;}else{mAaW++;}
					if(aaMatch(eaa, raa)){eAaC++;}else{eAaW++;}
					if(rac!=null){
						if(rac.equals(mac)){mAcC++;}else{mAcW++;}
						if(rac.equals(eac)){eAcC++;}else{eAcW++;}
					}
				}
			}
		}

		final long fn=refCount-tp;//oracle convention: unconsumed refs (may be <0 in pathological cases)
		outstream.println();
		outstream.println("=== "+type+" overlap-based metrics (reciprocal >50% overlap) ===");
		outstream.println("Overlap: \t"+tp+"tp\t"+fp+"fp\t"+fn+"fn\tprec="
			+Tools.format("%.1f%%", tp*100.0/Math.max(1, tp+fp))+"\trecall="
			+Tools.format("%.1f%%", tp*100.0/Math.max(1, refCount)));
		outstream.println("Boundary: \tboth-exact "+bothExact+"/"+tp+" ("
			+Tools.format("%.1f%%", bothExact*100.0/Math.max(1, tp))+"); 5p-exact "
			+Tools.format("%.1f%%", startExact*100.0/Math.max(1, tp))+", 3p-exact "
			+Tools.format("%.1f%%", stopExact*100.0/Math.max(1, tp)));
		outstream.println("Offsets: \t5p mean "+Tools.format("%+.2f", sumDs*1.0/Math.max(1, tp))
			+" |"+Tools.format("%.2f", absDs*1.0/Math.max(1, tp))+"|, 3p mean "
			+Tools.format("%+.2f", sumDe*1.0/Math.max(1, tp))+" |"+Tools.format("%.2f", absDe*1.0/Math.max(1, tp))+"|");
		if(isTrna){
			outstream.println("Extraction coverage: \t"+extracted+"/"+qryCount+" calls ("
				+Tools.format("%.1f%%", extracted*100.0/Math.max(1, qryCount))+")");
			outstream.println("Model-only AA: \t"+mAaC+"/"+(mAaC+mAaW)+" ("
				+Tools.format("%.2f%%", mAaC*100.0/Math.max(1, mAaC+mAaW))+")\tanticodon: "
				+mAcC+"/"+(mAcC+mAcW)+" ("+Tools.format("%.2f%%", mAcC*100.0/Math.max(1, mAcC+mAcW))+")");
			outstream.println("Extraction AA: \t"+eAaC+"/"+(eAaC+eAaW)+" ("
				+Tools.format("%.2f%%", eAaC*100.0/Math.max(1, eAaC+eAaW))+")\tanticodon: "
				+eAcC+"/"+(eAcC+eAcW)+" ("+Tools.format("%.2f%%", eAcC*100.0/Math.max(1, eAcC+eAcW))+")");
			outstream.println("Skipped (UNK ref): \t"+skipUnk);
		}
	}

	/** Amino-acid match with fMet credit (a Met call satisfies an fMet reference). */
	private static boolean aaMatch(String called, String ref){
		return called.equals(ref) || ("Met".equals(called) && "fMet".equals(ref));
	}

	/** Parses reference amino acid and (optional) anticodon triplet from a tRNA attribute string.
	 * Prefers the Note=tRNA-AA(ANTICODON) form (carries the triplet), else product=tRNA-AA.
	 * @return [aminoAcid, anticodonTripletOrNull]; aminoAcid is "UNK" if unrecognized */
	private static String[] refAmino(String attr){
		if(attr!=null){
			Matcher m=REF_NOTE.matcher(attr);
			if(m.find()){return new String[]{m.group(1), m.group(2)};}
			m=REF_PRODUCT.matcher(attr);
			if(m.find()){return new String[]{m.group(1), null};}
		}
		return new String[]{"UNK", null};
	}

	/** Parses the called anticodon/amino acid two ways: (a) model-name annotation
	 * (model:tRNA_consensus_XXX_c...) and (b) structural extraction (anticodon:XXX, DNA-&gt;RNA),
	 * falling back to the model name when no extraction attribute is present.
	 * @return [modelAA, modelAC, extractAA, extractAC, hasExtract("1"/"0")] */
	private static String[] calledAmino(String attr){
		String mac="NONE";
		if(attr!=null){
			Matcher m=CALL_MODEL.matcher(attr);
			if(m.find()){mac=m.group(1);}
		}
		final String maa=AC.getOrDefault(mac, mac);
		String eac, eaa, hasExtract;
		final Matcher e=(attr==null ? null : CALL_ANTI.matcher(attr));
		if(e!=null && e.find()){
			eac=e.group(1).replace('T', 'U');//DNA alphabet -> RNA
			eaa=AC.getOrDefault(eac, eac);
			hasExtract="1";
		}else{
			eac=mac; eaa=maa; hasExtract="0";
		}
		return new String[]{maa, mac, eaa, eac, hasExtract};
	}

	/** Builds the anticodon(RNA triplet)-&gt;amino-acid table verbatim from eval_anticodon3.sh
	 * (Aug 10 2026 revision). First writer wins on any duplicate triplet. */
	private static HashMap<String, String> buildAC(){
		HashMap<String, String> ac=new HashMap<String, String>();
		putAC(ac, "Phe", "AAA","GAA");
		putAC(ac, "Leu", "AAG","GAG","UAA","UAG","CAA","CAG");
		putAC(ac, "Ile", "AAU","GAU","UAU");
		putAC(ac, "Met", "CAU");
		putAC(ac, "Val", "AAC","GAC","UAC","CAC");
		putAC(ac, "Ser", "AGA","GGA","UGA","CGA","GCU","ACU");
		putAC(ac, "Pro", "AGG","GGG","UGG","CGG");
		putAC(ac, "Thr", "AGU","GGU","UGU","CGU");
		putAC(ac, "Ala", "AGC","GGC","UGC","CGC");
		putAC(ac, "Tyr", "AUA","GUA");
		putAC(ac, "His", "AUG","GUG");
		putAC(ac, "Gln", "UUG","CUG");
		putAC(ac, "Asn", "AUU","GUU");
		putAC(ac, "Lys", "UUU","CUU");
		putAC(ac, "Asp", "AUC","GUC");
		putAC(ac, "Glu", "UUC","CUC");
		putAC(ac, "Cys", "ACA","GCA");
		putAC(ac, "Trp", "CCA");
		putAC(ac, "Arg", "ACG","GCG","UCG","CCG","UCU","CCU");
		putAC(ac, "Gly", "ACC","GCC","UCC","CCC");
		putAC(ac, "Sec", "UCA");
		return ac;
	}

	/** Adds each anticodon triplet under an amino acid, keeping the first assignment. */
	private static void putAC(HashMap<String, String> map, String aa, String... acs){
		for(String a : acs){if(!map.containsKey(a)){map.put(a, aa);}}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in=null;
	private String ref=null;
	
	
	/*--------------------------------------------------------------*/

	private HashMap<StringNum, GffLine> lineMap;
	private HashMap<StringNum, Integer> startCountMap;
	private HashMap<StringNum, Integer> stopCountMap;
	
//	private HashMap<Integer, ArrayList<GffLine>> map;
//	private HashSet<Integer> stopSet;
//	private HashSet<Integer> startSet;
//	private HashSet<Integer> stopSetM;
//	private HashSet<Integer> startSetM;
	
	private long linesProcessed=0;
	private long linesOut=0;
	private long bytesProcessed=0;
	private long bytesOut=0;
	
	private long maxLines=Long.MAX_VALUE;

	private long falsePositiveStart=0;
	private long falsePositiveStop=0;
	private long truePositiveStart=0;
	private long truePositiveStop=0;
	private long falseNegativeStart=0;
	private long falseNegativeStop=0;
	
	private long falsePositiveStart2=0;
	private long falsePositiveStop2=0;
	private long truePositiveStart2=0;
	private long truePositiveStop2=0;
	
	private long refCount=0;
	private long queryCount=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private final FileFormat ffin;
	private final FileFormat ffref;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;

	/** If true, print the additional overlap-based ncRNA (tRNA/rRNA) metrics after the exact-match output. */
	private boolean ncrna=true;

	/** Anticodon(RNA triplet)->amino-acid table (see {@link #buildAC}). */
	private static final HashMap<String, String> AC=buildAC();

	/** Reference tRNA amino acid + anticodon triplet: Note=tRNA-AA(ANTICODON). */
	private static final Pattern REF_NOTE=Pattern.compile("Note=tRNA-(\\w+)\\((\\w+)\\)");
	/** Reference tRNA amino acid only: product=tRNA-AA. */
	private static final Pattern REF_PRODUCT=Pattern.compile("product=tRNA-(\\w+)");
	/** Called tRNA model-name anticodon: model:tRNA_consensus_XXX_c. */
	private static final Pattern CALL_MODEL=Pattern.compile("model:tRNA_consensus_(\\w+)_c");
	/** Called tRNA structural anticodon (DNA alphabet): anticodon:XXX. */
	private static final Pattern CALL_ANTI=Pattern.compile("anticodon:([ACGT]+)");

	/** Orders GffLines by ascending start (col4) for per-seqid overlap scans. */
	private static final java.util.Comparator<GffLine> START_COMP=new java.util.Comparator<GffLine>(){
		@Override
		public int compare(GffLine x, GffLine y){return Integer.compare(x.start, y.start);}
	};

}
