package jgi;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import shared.Timer;
import shared.Tools;
import stream.SamLine;

public class SplitSamFile {
	
	
	public static void main(String[] args){
		
		Timer t=new Timer();
		
		String in=args[0];
		String outF=args.length>1 ? args[1] : null;
		String outR=args.length>2 ? args[2] : null;
		String outU=args.length>3 ? args[3] : null;
		if(args.length>4){
			if(args[4].equalsIgnoreCase("header")){includeHeader=true;}
		}
		
		ByteFile tf=ByteFile.makeByteFile(in, false);
		
		Tools.testForDuplicateFiles(true, in, outF, outR, outU);
		Tools.testOutputFiles(true, false, false, outF, outR, outU);
		
		final ByteStreamWriter fStream, rStream, uStream;
		
		fStream=(outF==null ? null : new ByteStreamWriter(outF, true, false, true));
		rStream=(outR==null ? null : new ByteStreamWriter(outR, true, false, true));
		uStream=(outU==null ? null : new ByteStreamWriter(outU, true, false, true));

		if(fStream!=null){fStream.start();}
		if(rStream!=null){rStream.start();}
		if(uStream!=null){uStream.start();}
		
		long plus=0;
		long minus=0;
		long other=0;
		
		byte[] s=null;
		for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s.length>0){
				byte c=s[0];
				if(c=='@'){
					if(includeHeader){
						if(fStream!=null){fStream.println(s);}
						if(rStream!=null){rStream.println(s);}
						if(uStream!=null){uStream.println(s);}
					}
				}else{
					int flag=SamLine.parseFlagOnly(s);
					if(SamLine.mapped(flag)){
						if(SamLine.strand(flag)==0){
							if(fStream!=null){fStream.println(s);}
							plus++;
						}else{
							if(rStream!=null){rStream.println(s);}
							minus++;
						}
					}else{
						if(uStream!=null){uStream.println(s);}
						other++;
					}
				}
			}
		}
		tf.close();
		if(fStream!=null){fStream.poisonAndWait();}
		if(rStream!=null){rStream.poisonAndWait();}
		if(uStream!=null){uStream.poisonAndWait();}
		
		System.err.println("Total reads:   \t"+(plus+minus+other));
		System.err.println("Plus reads:    \t"+(plus));
		System.err.println("Minus reads:   \t"+(minus));
		System.err.println("Unmapped reads:\t"+(other));
		
		t.stop();
		
		System.err.println("Time:          \t"+t);
		
	}
	
	private static boolean includeHeader=false;
}
