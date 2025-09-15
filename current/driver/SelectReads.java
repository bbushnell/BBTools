package driver;

import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Shared;
import shared.Tools;
import stream.SamLine;

/**
 * 
 * Selects only reads with long deletions
 * 
 * @author Brian Bushnell
 * @date Jun 21, 2013
 *
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
		if(symbol=='='){symbol='M';}
		if(symbol=='X'){symbol='S';}
		if(symbol=='N'){symbol='D';}
		if(symbol=='S' || symbol=='H' || symbol=='P'){symbol='C';}
		
		final int index=Tools.indexOf(new char[] {'M','S','D','I','C'}, symbol);
		assert(index>=0) : "Symbol (3rd argument) must be M, S, D, I, C (for match string symbols) or M, =, X, D, N, I, S, H, P (for cigar symbols).";
		
		TextFile tf=new TextFile(args[0], true);
		TextStreamWriter tsw=new TextStreamWriter(args[1], false, false, true);
		tsw.start();
		
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.charAt(0)=='@'){
				tsw.println(line);
			}else{
				if((reads=reads-1)<0){break;}
				SamLine sl=new SamLine(line);
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
