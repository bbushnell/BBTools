package driver;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import shared.Shared;
import shared.Tools;
import stream.SamLine;

/**
 * Selects reads from SAM files containing specific CIGAR operations above a minimum length.
 * Filters reads based on match, substitution, deletion, insertion, or clipping operations.
 * Processes SAM format input and outputs matching reads in SAM format.
 *
 * @author Brian Bushnell
 * @date Jun 21, 2013
 */
public final class SelectReads {
	
	public static void main(String[] args){
		
		assert(args.length>=2) : "Need 2 file names: <input> <output>";
		assert(!args[0].equalsIgnoreCase(args[1])) : "File names must be different.";
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		
		int minlen=1;
		long reads=Long.MAX_VALUE;
		char symbol='D';
		if(args.length>2){symbol=(char)args[2].charAt(0);}
		if(args.length>3){minlen=Integer.parseInt(args[3]);}
		if(args.length>4){reads=Parse.parseKMG(args[4]);}
		
		symbol=Tools.toUpperCase(symbol);
		//TODO: Possible bug [driver/SelectReads#001] LOW/dev (no .sh, no callers): these are SEQUENTIAL ifs, not else-if, so
		//they CASCADE. cigar 'X' → 'S' (L40) then immediately 'S' → 'C' (L42), and a direct 'S' also → 'C'. The target index
		//array below is {'M','S','D','I','C'}, where index 1 = 'S' (substitution). But NO input can ever produce a final 'S':
		//both 'X' and 'S' collapse to 'C'. → substitution selection is UNREACHABLE; msdic[1] is never queryable. Intent is
		//ambiguous (match-string 'S'=substitution vs cigar 'S'=softclip), so ESCALATE rather than blind-fix (else-if chain, or
		//pick which 'S' wins). Dead dev tool → LOW.
		if(symbol=='='){symbol='M';}
		if(symbol=='X'){symbol='S';}
		if(symbol=='N'){symbol='D';}
		if(symbol=='S' || symbol=='H' || symbol=='P'){symbol='C';}

		final int index=Tools.indexOf(new char[] {'M','S','D','I','C'}, symbol);
		assert(index>=0) : "Symbol (3rd argument) must be M, S, D, I, C (for match string symbols) or M, =, X, D, N, I, S, H, P (for cigar symbols).";
		
		ByteFile tf=ByteFile.makeByteFile(args[0], true);
		LineParser1 lp=new LineParser1('\t');
		ByteStreamWriter tsw=new ByteStreamWriter(args[1], false, false, true);
		tsw.start();
		
		for(byte[] line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line[0]=='@'){
				tsw.println(line);
			}else{
				if((reads=reads-1)<0){break;}
				SamLine sl=new SamLine(lp.set(line));
				if(testLine(sl, minlen, index)){
					tsw.println(line);
				}
			}
		}
		tf.close();
		tsw.poisonAndWait();
		
	}
	
	
	private static boolean testLine(SamLine sl, int minlen, int index){
		assert(sl!=null);
		if(!sl.mapped() || sl.cigar==null){return false;}
		int[] msdic=sl.cigarToMdsiMax(sl.cigar);
		return (msdic!=null && msdic[index]>=minlen);
	}
	
}
