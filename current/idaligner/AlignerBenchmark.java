package idaligner;

import java.util.ArrayList;
import java.util.Arrays;
import align2.MultiStateAligner11ts;
import dna.AminoAcid;
import fileIO.TextFile;
import structures.ByteBuilder;
import structures.LongList;

/**
 * Test harness for comparing aligners against ground-truth synthetic reads.
 * Generates reads with randomreads.sh, parses SYN headers for truth coordinates,
 * extracts padded reference windows, and aligns with any IDAligner or the 4SA (MSA11ts).
 *
 * Usage:
 *   java idaligner.AlignerBenchmark ref=phix.fa in=reads.fq [pad=20]
 *       [aligners=ScrabbleAligner,ScrabbleAffine,MSA11ts] [tracescrabble=t]
 *
 * @author Ady
 * @date July 1, 2026
 */
public class AlignerBenchmark{

	public static void main(String[] args){
		String refPath=null;
		String readsPath=null;
		int padding=20;
		String[] alignerNames={"ScrabbleAligner", "ScrabbleAffine", "MSA11ts"};
		boolean traceScrabble=true;
		boolean parseCigar=true;

		for(String arg : args){
			String[] kv=arg.split("=", 2);
			String key=kv[0].toLowerCase();
			String val=kv.length>1 ? kv[1] : "";
			if(key.equals("ref")){refPath=val;}
			else if(key.equals("in") || key.equals("reads")){readsPath=val;}
			else if(key.equals("pad") || key.equals("padding")){padding=Integer.parseInt(val);}
			else if(key.equals("aligners")){alignerNames=val.split(",");}
			else if(key.equals("tracescrabble") || key.equals("trace")){traceScrabble=val.startsWith("t");}
			else if(key.equals("parsecigar") || key.equals("cigar")){parseCigar=parse.Parse.parseBoolean(val);}
		}

		if(refPath==null || readsPath==null){
			System.err.println("Usage: AlignerBenchmark ref=<fasta> in=<fastq> [pad=20]");
			System.err.println("       [aligners=ScrabbleAligner,ScrabbleAffine,MSA11ts]");
			System.err.println("       [trace=t]");
			System.err.println();
			System.err.println("Generate reads first:");
			System.err.println("  randomreads.sh ref=phix.fa out=reads.fq reads=100 len=150 \\");
			System.err.println("    snprate=0.01 insrate=0.005 delrate=0.005 adderrors=f");
			System.exit(1);
		}

		byte[] reference=loadReference(refPath);
		System.err.println("Loaded reference: "+reference.length+" bp");

		// Build aligners: each entry is either an IDAligner or the string "MSA11ts"
		ArrayList<Object> alignerList=new ArrayList<>();
		ArrayList<String> nameList=new ArrayList<>();
		for(String name : alignerNames){
			if(name.equalsIgnoreCase("MSA11ts") || name.equalsIgnoreCase("MSA") || name.equalsIgnoreCase("4SA")){
				alignerList.add("MSA11ts");
				nameList.add("MSA11ts");
			}else{
				IDAligner a=makeAligner(name);
				if(a!=null){alignerList.add(a); nameList.add(a.name());}
				else{System.err.println("Unknown aligner: "+name);}
			}
		}

		// Print header
		ByteBuilder header=new ByteBuilder();
		header.append("readID\tstrand\ttrueStart\ttrueStop\treadLen");
		if(parseCigar){header.append("\ttruthCigar");}
		for(int i=0; i<nameList.size(); i++){
			String n=nameList.get(i);
			header.append('\t').append(n).append("_score");
			header.append('\t').append(n).append("_start");
			header.append('\t').append(n).append("_stop");
			header.append('\t').append(n).append("_match");
			if(parseCigar){header.append('\t').append(n).append("_cigarMatch");}
		}
		System.out.println(header);

		// Per-aligner stats for CIGAR comparison
		int[][] cigarStats=new int[nameList.size()][4]; // [exact, opsMatch, opsMismatch, noCompare]

		// Process reads
		TextFile tf=new TextFile(readsPath);
		String line;
		int readCount=0;
		int cigarCount=0;
		LongList trace=traceScrabble ? new LongList(8192) : null;
		ByteBuilder bb=new ByteBuilder(512);
		MultiStateAligner11ts msa=null;

		while((line=tf.nextLine())!=null){
			if(!line.startsWith("@SYN_") && !line.startsWith("@syn_")){continue;}
			String headerStr=line.substring(1); // remove @
			String seqLine=tf.nextLine();
			tf.nextLine(); // + line
			tf.nextLine(); // quality line
			if(seqLine==null){break;}

			byte[] query=seqLine.trim().getBytes();
			SynHeader sh=parseHeader(headerStr);

			// Extract padded reference window
			int winStart=Math.max(0, sh.start-padding);
			int winStop=Math.min(reference.length-1, sh.stop+padding);
			byte[] refWindow=Arrays.copyOfRange(reference, winStart, winStop+1);

			// Reverse complement query if minus strand
			byte[] alignQuery=query;
			if(sh.strand==1){
				alignQuery=AminoAcid.reverseComplementBases(query);
			}

			// Expand truth CIGAR to per-base match string if present
			byte[] truthExpanded=null;
			if(parseCigar && sh.truthMatch!=null){
				truthExpanded=cigarToMatchString(sh.truthMatch);
				cigarCount++;
			}

			ByteBuilder row=new ByteBuilder();
			row.append(sh.id).append('\t').append(sh.strand==0 ? '+' : '-');
			row.append('\t').append(sh.start).append('\t').append(sh.stop);
			row.append('\t').append(query.length);
			if(parseCigar){
				row.append('\t').append(sh.truthMatch==null ? "." : new String(sh.truthMatch));
			}

			int[] pos=new int[3];

			for(int i=0; i<nameList.size(); i++){
				String name=nameList.get(i);

				byte[] alignerMatch=null;
				int alignerStart=-1, alignerStop=-1;

				if(alignerList.get(i) instanceof String){
					// MSA11ts alignment: fillUnlimited returns {rows, maxCol, maxState, maxScore}
					if(msa==null || msa.maxRows<alignQuery.length || msa.maxColumns<refWindow.length){
						msa=new MultiStateAligner11ts(alignQuery.length, refWindow.length+10);
					}
					int[] max=msa.fillUnlimited(alignQuery, refWindow, 0, refWindow.length-1, null);
					int[] scoreArr=msa.score(alignQuery, refWindow, 0, refWindow.length-1,
						max[0], max[1], max[2], false);
					alignerMatch=msa.traceback(alignQuery, refWindow, 0, refWindow.length-1,
						max[0], max[1], max[2], false);
					alignerStart=winStart+scoreArr[1];
					alignerStop=winStart+scoreArr[2];
					row.append('\t').append(scoreArr[0]);
					row.append('\t').append(alignerStart);
					row.append('\t').append(alignerStop);
					row.append('\t').append(alignerMatch==null ? "." : new String(alignerMatch));

				}else if(name.equals("Scrabble")){
					AlignmentStats st=new AlignmentStats(true);
					ScrabbleAligner.alignAndTraceStatic(alignQuery, refWindow, st);
					alignerMatch=st.matchString;
					alignerStart=winStart+st.rStart;
					alignerStop=winStart+st.rStop;
					row.append('\t').append(st.score);
					row.append('\t').append(alignerStart);
					row.append('\t').append(alignerStop);
					row.append('\t').append(alignerMatch==null ? "." : new String(alignerMatch));

				}else if(name.equals("ScrabbleAffine") && trace!=null){
					IDAligner aligner=(IDAligner)alignerList.get(i);
					aligner.align(alignQuery, refWindow, pos);
					alignerStart=winStart+pos[0];
					alignerStop=winStart+pos[1];
					row.append('\t').append(pos[2]).append('\t').append(alignerStart).append('\t').append(alignerStop);
					ScrabbleAffine.alignWithTrace(alignQuery, refWindow, pos, trace);
					alignerMatch=TracerAffine.traceback(trace, alignQuery, refWindow,
						alignQuery.length, pos[1]+1, bb);
					row.append('\t').append(new String(alignerMatch));

				}else{
					IDAligner aligner=(IDAligner)alignerList.get(i);
					aligner.align(alignQuery, refWindow, pos);
					alignerStart=winStart+pos[0];
					alignerStop=winStart+pos[1];
					row.append('\t').append(pos[2]).append('\t').append(alignerStart).append('\t').append(alignerStop);
					row.append('\t').append('.');
				}

				// Compare against truth CIGAR if available
				if(parseCigar && truthExpanded!=null && alignerMatch!=null){
					String result=compareToCigar(truthExpanded, alignerMatch,
						sh.start, sh.stop, alignerStart, alignerStop);
					row.append('\t').append(result);
					if(result.equals("EXACT")){cigarStats[i][0]++;}
					else if(result.startsWith("OPS_OK")){cigarStats[i][1]++;}
					else if(result.startsWith("DIFF")){cigarStats[i][2]++;}
					else{cigarStats[i][3]++;}
				}else if(parseCigar){
					row.append('\t').append('.');
				}
			}

			System.out.println(row);
			readCount++;
		}
		tf.close();
		System.err.println("Processed "+readCount+" reads.");

		// Print CIGAR comparison summary
		if(parseCigar && cigarCount>0){
			System.err.println("\n--- CIGAR comparison ("+cigarCount+" reads with truth CIGARs) ---");
			for(int i=0; i<nameList.size(); i++){
				int exact=cigarStats[i][0], opsOk=cigarStats[i][1];
				int diff=cigarStats[i][2], noCmp=cigarStats[i][3];
				int compared=exact+opsOk+diff;
				if(compared==0){continue;}
				System.err.println(nameList.get(i)+": exact="+exact+"/"+compared
					+" opsMatch="+opsOk+" diff="+diff+" noCompare="+noCmp);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Reference loading      ----------------*/
	/*--------------------------------------------------------------*/

	static byte[] loadReference(String path){
		TextFile tf=new TextFile(path);
		ByteBuilder bb=new ByteBuilder(8192);
		String line;
		while((line=tf.nextLine())!=null){
			if(line.startsWith(">")){continue;}
			bb.append(line.trim().toUpperCase());
		}
		tf.close();
		return bb.toBytes();
	}

	/*--------------------------------------------------------------*/
	/*----------------      Header parsing         ----------------*/
	/*--------------------------------------------------------------*/

	/** Parse a SYN_ header into its components. */
	static SynHeader parseHeader(String header){
		// Format: SYN_id_start_stop_insert_strand_bbstart_bbchrom_match_rname [pairnum:]
		// Remove trailing " 1:" or "/1" pairnum suffix
		int space=header.indexOf(' ');
		if(space>0){header=header.substring(0, space);}
		// Remove "SYN_" prefix
		String body=header.substring(4);
		String[] parts=body.split("_", 10);
		SynHeader sh=new SynHeader();
		sh.id=Long.parseLong(parts[0]);
		sh.start=Integer.parseInt(parts[1]);
		sh.stop=Integer.parseInt(parts[2]);
		sh.strand=parts[4].equals("+") ? 0 : 1;
		if(parts.length>7 && !parts[7].equals(".")){
			sh.truthMatch=parts[7].getBytes();
		}
		return sh;
	}

	static class SynHeader{
		long id;
		int start, stop, strand;
		byte[] truthMatch;
	}

	/*--------------------------------------------------------------*/
	/*----------------      CIGAR comparison       ----------------*/
	/*--------------------------------------------------------------*/

	/** Expand run-length CIGAR (70=5I20=) to per-base match string (mmm...IIIIImmm...). */
	static byte[] cigarToMatchString(byte[] cigar){
		ByteBuilder bb=new ByteBuilder(cigar.length*2);
		int num=0;
		for(int i=0; i<cigar.length; i++){
			byte c=cigar[i];
			if(c>='0' && c<='9'){
				num=num*10+(c-'0');
			}else{
				byte op;
				if(c=='='){op='m';}
				else if(c=='X'){op='S';}
				else if(c=='I'){op='I';}
				else if(c=='D'){op='D';}
				else{op='N';}
				for(int j=0; j<num; j++){bb.append(op);}
				num=0;
			}
		}
		return bb.toBytes();
	}

	/** Count operations in a per-base match string. Returns {matches, subs, ins, dels}. */
	static int[] countOps(byte[] match){
		int m=0, s=0, ins=0, del=0;
		for(byte b : match){
			if(b=='m'){m++;}
			else if(b=='S'){s++;}
			else if(b=='I'){ins++;}
			else if(b=='D'){del++;}
		}
		return new int[]{m, s, ins, del};
	}

	/** Count indel run opens in a per-base match string. Returns {delOpens, insOpens}. */
	static int[] countOpens(byte[] match){
		int dOpens=0, iOpens=0;
		for(int i=0; i<match.length; i++){
			if(match[i]=='D' && (i==0 || match[i-1]!='D')){dOpens++;}
			if(match[i]=='I' && (i==0 || match[i-1]!='I')){iOpens++;}
		}
		return new int[]{dOpens, iOpens};
	}

	/**
	 * Compare an aligner's match string to the truth CIGAR (both in per-base format).
	 * Returns: EXACT, OPS_OK(counts match but positions differ),
	 * DIFF(m:+2/-1,I:+0/-0,D:+0/-1), or POS_MISMATCH.
	 */
	static String compareToCigar(byte[] truth, byte[] aligner,
			int trueStart, int trueStop, int alignerStart, int alignerStop){
		if(trueStart!=alignerStart || trueStop!=alignerStop){
			return "POS_MISMATCH";
		}
		if(Arrays.equals(truth, aligner)){return "EXACT";}

		int[] tOps=countOps(truth), aOps=countOps(aligner);
		int[] tOpens=countOpens(truth), aOpens=countOpens(aligner);
		boolean opsMatch=(tOps[0]==aOps[0] && tOps[1]==aOps[1]
			&& tOps[2]==aOps[2] && tOps[3]==aOps[3]);
		if(opsMatch){return "OPS_OK(opens:D="+aOpens[0]+"v"+tOpens[0]+",I="+aOpens[1]+"v"+tOpens[1]+")";}

		return "DIFF(m:"+delta(tOps[0],aOps[0])+",S:"+delta(tOps[1],aOps[1])
			+",I:"+delta(tOps[2],aOps[2])+",D:"+delta(tOps[3],aOps[3])
			+",opens:D="+aOpens[0]+"v"+tOpens[0]+",I="+aOpens[1]+"v"+tOpens[1]+")";
	}

	private static String delta(int truth, int aligner){
		int d=aligner-truth;
		return d==0 ? "0" : (d>0 ? "+"+d : ""+d);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Aligner factory        ----------------*/
	/*--------------------------------------------------------------*/

	static IDAligner makeAligner(String name){
		name=name.toLowerCase();
		if(name.equals("scrabblealigner") || name.equals("scrabble")){return new ScrabbleAligner();}
		if(name.equals("scrabbleaffine") || name.equals("affine")){return new ScrabbleAffine();}
		if(name.equals("glocal+5") || name.equals("glocal5") || name.equals("glocalplus5")){return new GlocalPlusAligner5();}
		if(name.equals("glocal+4") || name.equals("glocal4") || name.equals("glocalplus4")){return new GlocalPlusAligner4();}
		if(name.equals("glocal+3") || name.equals("glocal3") || name.equals("glocalplus3")){return new GlocalPlusAligner3();}
		if(name.equals("glocal+2") || name.equals("glocal2") || name.equals("glocalplus2")){return new GlocalPlusAligner2();}
		if(name.equals("glocal+") || name.equals("glocal1") || name.equals("glocalplus")){return new GlocalPlusAligner();}
		if(name.equals("glocal") || name.equals("glocalaligner")){return new GlocalAligner();}
		if(name.equals("glocalint")){return new GlocalAlignerInt();}
		if(name.equals("glocalconcise")){return new GlocalAlignerConcise();}
		if(name.equals("scrabble2") || name.equals("scrabblealigner2")){return new ScrabbleAligner2();}
		return null;
	}

}
