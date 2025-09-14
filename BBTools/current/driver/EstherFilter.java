package driver;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Jul 29, 2014
 *
 */
public class EstherFilter {
	
	public static void main(String[] args){
		String query=args[0];
		String ref=args[1];
		float cutoff=Float.parseFloat(args[2]);
		boolean outputFasta=false;
		if(args.length>3 && args[3].equalsIgnoreCase("fasta")){
			outputFasta=true;
		}
		String command="blastall -p blastn -i "+query+" -d "+ref+" -e 0.00001 -m 8";
		
		ReadWrite.FORCE_KILL=true;

//		InputStream is=ReadWrite.getInputStreamFromProcess("stdin", command, false);
//		InputStream is=ReadWrite.getInputStreamFromProcess("", command, false);
		InputStream is=ReadWrite.getInputStreamFromProcess("foo", command, false, false, true);
		
		InputStreamReader isr=new InputStreamReader(is);
		BufferedReader b=new BufferedReader(isr, 32768);
		
//		System.out.println("Before");
		
		if(outputFasta){
//			System.out.println("A");
			processToFasta(b, cutoff, query);
		}else{
//			System.out.println("B");
			processToNames(b, cutoff);
		}
		
//		System.out.println("Finished");

//		ReadWrite.finishReading(is, "stdin", true, b, isr);
//		ReadWrite.finishReading(is, "", true, b, isr);
		ReadWrite.finishReading(is, "foo", true, b, isr);
		
	}
	
	public static void processToFasta(BufferedReader b, float cutoff, String query){
		String s=null;
		
		ArrayList<String> names=new ArrayList<String>();
//		System.out.println("Reading line 0");
		try {
			s=b.readLine();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
//		System.out.println("Starting");
		String prev="";
		
		while(s!=null){
			String[] split=s.split("\t");
			float value=0;
			try {
				value=Float.parseFloat(split[11].trim());
			} catch (NumberFormatException e) {
				e.printStackTrace();
//				System.err.println("Bad line:\n"+s);
			}
			if(value>=cutoff){
				if(!prev.equals(split[0])){
					prev=split[0];
					names.add(split[0]);
				}
			}
//			System.out.println("Reading line");
			try {
				s=b.readLine();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		outputFasta(query, names);
	}
	
	public static void processToNames(BufferedReader b, float cutoff){
		String s=null;
//		System.out.println("Reading line 0");
		try {
			s=b.readLine();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
//		System.out.println("Starting");
		String prev="";
		while(s!=null){
			String[] split=s.split("\t");
			float value=0;
			try {
				value=Float.parseFloat(split[11].trim());
			} catch (NumberFormatException e) {
				e.printStackTrace();
//				System.err.println("Bad line:\n"+s);
			}
			if(value>=cutoff){
				System.out.println(split[0]);
			}
//			System.out.println("Reading line");
			try {
				s=b.readLine();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public static void outputFasta(String fname, ArrayList<String> names){
		
		Shared.sort(names);
		
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, false, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1L, false, ff, null);
		cris.start(); //4567
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		/* Iterate through read lists from the input stream */
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			
			for(Read r : reads){
				if(Collections.binarySearch(names, r.id)>=0){
					System.out.println(r.toFasta(70));
				}
			}
			
			/* Dispose of the old list and fetch a new one */
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		/* Cleanup */
		cris.returnList(ln);
	}
	
	
}
