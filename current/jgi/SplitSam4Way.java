package jgi;

import java.io.PrintStream;

import fileIO.FileFormat;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.SamLine;

/**
 * @author Brian Bushnell
 * @date Jul 23, 2013
 *
 */
public class SplitSam4Way {
	
	public static void main(String[] args){
		SplitSam4Way x=new SplitSam4Way(args);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	private void printOptions(){
		outstream.println("Syntax:\n");
		outstream.println("java -ea -Xmx128m -cp <path> jgi.SplitSam4Way <input> <out plus> <out minus> <out chimeric> <out unmapped>");
		outstream.println("If you do not want one of the output files, use the word 'null'.\n");
	}
	
	public SplitSam4Way(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		if(args==null || args.length!=5){
			printOptions();
			System.exit(1);
		}
		
		Timer t=new Timer();
		long reads=0, bases=0;
		long preads=0, mreads=0, creads=0, ureads=0;
		
		String fin=args[0];
		String fplus=args[1];
		String fminus=args[2];
		String fchimeric=args[3];
		String funmapped=args[4];
		
		TextFile tf=new TextFile(fin, true);
		TextStreamWriter plus=("null".equalsIgnoreCase(fplus) ? null : new TextStreamWriter(fplus, true, false, true, FileFormat.SAM));
		TextStreamWriter minus=("null".equalsIgnoreCase(fminus) ? null : new TextStreamWriter(fminus, true, false, true, FileFormat.SAM));
		TextStreamWriter chimeric=("null".equalsIgnoreCase(fchimeric) ? null : new TextStreamWriter(fchimeric, true, false, true, FileFormat.SAM));
		TextStreamWriter unmapped=("null".equalsIgnoreCase(funmapped) ? null : new TextStreamWriter(funmapped, true, false, true, FileFormat.SAM));

		if(plus!=null){plus.start();}
		if(minus!=null){minus.start();}
		if(chimeric!=null){chimeric.start();}
		if(unmapped!=null){unmapped.start();}
		
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.charAt(0)=='@'){
				if(plus!=null){plus.println(line);}
				if(minus!=null){minus.println(line);}
				if(chimeric!=null){chimeric.println(line);}
				if(unmapped!=null){unmapped.println(line);}
			}else{
				SamLine sl=new SamLine(line);
				reads++;
//				bases+=sl.seq.length();
				bases+=sl.seq.length;
				
				if(!sl.mapped() || !sl.nextMapped() || !sl.hasMate() || !sl.primary()){
					if(unmapped!=null){unmapped.println(line);}
					ureads++;
//					System.out.println("unmapped: "+sl.mapped()+", "+sl.nextMapped()+", "+sl.hasMate()+", "+!sl.primary());
				}else if(!sl.pairedOnSameChrom() || sl.strand()==sl.nextStrand()){
					if(chimeric!=null){chimeric.println(line);}
					creads++;
//					System.out.println("chimeric: "+sl.pairedOnSameChrom()+", "+(sl.strand()==sl.nextStrand())+", "+sl.strand()+", "+sl.nextStrand()+", "+new String(sl.rname())+", "+new String(sl.rnext()));
				}else if((sl.firstFragment() ? sl.strand() : sl.nextStrand())==Shared.PLUS){
					if(plus!=null){plus.println(line);}
					preads++;
				}else if((sl.firstFragment() ? sl.strand() : sl.nextStrand())==Shared.MINUS){
					if(minus!=null){minus.println(line);}
					mreads++;
				}else{
					throw new RuntimeException("Unhandled case: "+sl.firstFragment()+", "+sl.lastFragment()+", "+sl.strand()+", "+sl.nextStrand()+"\n"+sl+"\n");
				}
			}
		}
		
		if(plus!=null){plus.poisonAndWait();}
		if(minus!=null){minus.poisonAndWait();}
		if(chimeric!=null){chimeric.poisonAndWait();}
		if(unmapped!=null){unmapped.poisonAndWait();}
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
		outstream.println("Plus Reads:         "+preads);
		outstream.println("Minus Reads:        "+mreads);
		outstream.println("Chimeric Reads:     "+creads);
		outstream.println("Unmapped Reads:     "+ureads);
	}
	
	private PrintStream outstream=System.err;
	
}
