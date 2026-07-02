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

		for(String arg : args){
			String[] kv=arg.split("=", 2);
			String key=kv[0].toLowerCase();
			String val=kv.length>1 ? kv[1] : "";
			if(key.equals("ref")){refPath=val;}
			else if(key.equals("in") || key.equals("reads")){readsPath=val;}
			else if(key.equals("pad") || key.equals("padding")){padding=Integer.parseInt(val);}
			else if(key.equals("aligners")){alignerNames=val.split(",");}
			else if(key.equals("tracescrabble") || key.equals("trace")){traceScrabble=val.startsWith("t");}
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
		for(int i=0; i<nameList.size(); i++){
			String n=nameList.get(i);
			header.append('\t').append(n).append("_score");
			header.append('\t').append(n).append("_start");
			header.append('\t').append(n).append("_stop");
			header.append('\t').append(n).append("_match");
		}
		System.out.println(header);

		// Process reads
		TextFile tf=new TextFile(readsPath);
		String line;
		int readCount=0;
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

			ByteBuilder row=new ByteBuilder();
			row.append(sh.id).append('\t').append(sh.strand==0 ? '+' : '-');
			row.append('\t').append(sh.start).append('\t').append(sh.stop);
			row.append('\t').append(query.length);

			int[] pos=new int[3];

			for(int i=0; i<nameList.size(); i++){
				String name=nameList.get(i);

				if(alignerList.get(i) instanceof String){
					// MSA11ts alignment: fillUnlimited returns {rows, maxCol, maxState, maxScore}
					if(msa==null || msa.maxRows<alignQuery.length || msa.maxColumns<refWindow.length){
						msa=new MultiStateAligner11ts(alignQuery.length, refWindow.length+10);
					}
					int[] max=msa.fillUnlimited(alignQuery, refWindow, 0, refWindow.length-1, null);
					int[] scoreArr=msa.score(alignQuery, refWindow, 0, refWindow.length-1,
						max[0], max[1], max[2], false);
					byte[] msaMatch=msa.traceback(alignQuery, refWindow, 0, refWindow.length-1,
						max[0], max[1], max[2], false);

					row.append('\t').append(scoreArr[0]);
					row.append('\t').append(winStart+scoreArr[1]);
					row.append('\t').append(winStart+scoreArr[2]);
					row.append('\t').append(msaMatch==null ? "." : new String(msaMatch));

				}else if(name.equals("Scrabble")){
					// ScrabbleAligner with trace — uses alignAndTraceStatic for match strings
					AlignmentStats st=new AlignmentStats(true);
					ScrabbleAligner.alignAndTraceStatic(alignQuery, refWindow, st);
					row.append('\t').append(st.score);
					row.append('\t').append(winStart+st.rStart);
					row.append('\t').append(winStart+st.rStop);
					row.append('\t').append(st.matchString==null ? "." : new String(st.matchString));

				}else if(name.equals("ScrabbleAffine") && trace!=null){
					// ScrabbleAffine with TracerAffine for match strings
					IDAligner aligner=(IDAligner)alignerList.get(i);
					aligner.align(alignQuery, refWindow, pos);
					int absStart=winStart+pos[0];
					int absStop=winStart+pos[1];
					row.append('\t').append(pos[2]).append('\t').append(absStart).append('\t').append(absStop);
					ScrabbleAffine.alignWithTrace(alignQuery, refWindow, pos, trace);
					byte[] match=TracerAffine.traceback(trace, alignQuery, refWindow,
						alignQuery.length, pos[1]+1, bb);
					row.append('\t').append(new String(match));

				}else{
					// Generic IDAligner — no match string
					IDAligner aligner=(IDAligner)alignerList.get(i);
					aligner.align(alignQuery, refWindow, pos);
					int absStart=winStart+pos[0];
					int absStop=winStart+pos[1];
					row.append('\t').append(pos[2]).append('\t').append(absStart).append('\t').append(absStop);
					row.append('\t').append('.');
				}
			}

			System.out.println(row);
			readCount++;
		}
		tf.close();
		System.err.println("Processed "+readCount+" reads.");
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
