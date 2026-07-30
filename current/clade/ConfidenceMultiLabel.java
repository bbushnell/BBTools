package clade;

import java.io.PrintStream;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Converts a single-continuous-label confidence vector file (#dims 48 1: 48 raw features + one
 * continuous depth label in [0,1]) into a MULTI-BOOLEAN-label file (#dims 48 8: same 48 features +
 * 8 per-level booleans, "correct at level L"). This is the training/validation data for a SINGLE net
 * with 8 boolean outputs -- the middle option between one scalar net (bottlenecked) and eight separate
 * nets (expensive): one forward pass, 8 heads, so each level gets its own free parameter.
 *
 * "Correct at level L" = trueDepth >= encodeLevel(L), where encodeLevel(level)=(11-level)/9 (the same
 * ladder ConfidenceVectorizer/CladeConfidence use). Output order is species..domain (superkingdom
 * skipped) to match the shipped boolean bundle and clade.ConfidenceHeadToHead.
 *
 * Usage: confidencemultilabel.sh in=<vec_1.tsv> out=<vec_8.tsv> [maxrows=0] [stride=1]
 *   maxrows=  cap output rows (0 = all)
 *   stride=   keep 1 in N input rows (1 = all); use to subsample for an in-RAM trainer.
 *
 * @author Noire
 */
public class ConfidenceMultiLabel {

	static final int[] LCODE={2,3,4,5,6,7,8,10};   //species,genus,family,order,class,phylum,kingdom,domain
	static final int NL=LCODE.length;
	static final float EPS=1e-4f;

	public static void main(String[] args){
		String in=null, out=null; long maxrows=0; int stride=1, indims=48;
		for(String a : args){
			final String[] s=a.split("=", 2);
			final String k=s[0].toLowerCase(), v=(s.length>1 ? s[1] : null);
			if(k.equals("in")){in=v;}
			else if(k.equals("out")){out=v;}
			else if(k.equals("maxrows")){maxrows=Long.parseLong(v);}
			else if(k.equals("stride")){stride=Integer.parseInt(v);}
			else if(k.equals("indims")){indims=Integer.parseInt(v);}
			else{throw new RuntimeException("Unknown parameter: "+a);}
		}
		if(in==null || out==null){throw new RuntimeException("Required: in= out=");}

		final float[] thresh=new float[NL];
		for(int L=0; L<NL; L++){thresh[L]=(11-LCODE[L])/9f;}

		final ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(in, FileFormat.TEXT, null, true, false));
		final ByteStreamWriter bsw=new ByteStreamWriter(FileFormat.testOutput(out, FileFormat.TEXT, null, true, true, false, false));
		bsw.start();
		bsw.println(new ByteBuilder("#dims\t"+indims+"\t8"));

		final ByteBuilder bb=new ByteBuilder(256);
		byte[] line; long seen=0, written=0, bad=0;
		while((line=bf.nextLine())!=null){
			if(line.length==0 || line[0]=='#'){continue;}
			if((seen++ % stride)!=0){continue;}
			final String[] f=new String(line).split("\t");
			if(f.length!=indims+1){bad++; continue;}
			final float depth;
			try{depth=Float.parseFloat(f[indims]);}catch(Exception e){bad++; continue;}
			bb.clear();
			for(int i=0; i<indims; i++){bb.append(f[i]).tab();}      //indims features verbatim
			for(int L=0; L<NL; L++){
				bb.append(depth>=thresh[L]-EPS ? 1 : 0);
				if(L<NL-1){bb.tab();}
			}
			bsw.println(bb);
			written++;
			if(maxrows>0 && written>=maxrows){break;}
		}
		bf.close();
		bsw.poisonAndWait();
		final PrintStream e=System.err;
		e.println("in rows seen="+seen+"  written="+written+"  malformed="+bad);
	}
}
