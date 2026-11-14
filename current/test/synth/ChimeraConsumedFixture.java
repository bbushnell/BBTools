package test.synth;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;

/** Writes and checks MakeChimeras consumed-output benchmark FASTA fixtures.
 * @author Jean
 */
public final class ChimeraConsumedFixture {

	private ChimeraConsumedFixture(){}

	public static void main(String[] args) throws IOException{
		if(args.length<1){throw new IllegalArgumentException(usage());}
		if("write-inputs".equals(args[0])){
			if(args.length!=5){throw new IllegalArgumentException(usage());}
			writeInputs(new File(args[1]), Integer.parseInt(args[2]),
					Integer.parseInt(args[3]), Integer.parseInt(args[4]));
		}else if("check-fasta".equals(args[0])){
			if(args.length!=6){throw new IllegalArgumentException(usage());}
			checkFasta(args[1], args[2], new File(args[3]),
					Integer.parseInt(args[4]), Integer.parseInt(args[5]));
		}else{
			throw new IllegalArgumentException(usage());
		}
	}

	private static String usage(){
		return "Usage: ChimeraConsumedFixture write-inputs <dir> <count> <eligible length> <short length> "
				+ "| check-fasta <label> <workload> <path> <expected records> <expected length>";
	}

	private static void writeInputs(File dir, int count, int eligibleLength, int shortLength)
			throws IOException{
		if(!dir.isDirectory() && !dir.mkdirs()){
			throw new IOException("Could not create " + dir);
		}
		try(BufferedWriter out=writer(new File(dir, "all.fa"))){
			for(int i=1; i<=count; i++){
				writeRecord(out, "all", i, eligibleLength, "ACGT", 0);
			}
		}
		try(BufferedWriter out=writer(new File(dir, "mixed.fa"))){
			for(int i=1; i<=count; i++){
				writeRecord(out, "short", i, shortLength, "TGCA", 0);
			}
			for(int i=1; i<=count; i++){
				writeRecord(out, "eligible", i, eligibleLength, "TGCA", 1);
			}
		}
	}

	private static BufferedWriter writer(File file) throws IOException{
		return new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file),
				StandardCharsets.US_ASCII));
	}

	private static void writeRecord(BufferedWriter out, String id, int index, int length,
			String bases, int offset) throws IOException{
		out.write('>');
		out.write(id);
		out.write('_');
		writeFourDigits(out, index);
		out.write('\n');
		for(int i=0; i<length; i++){
			out.write(bases.charAt((offset+index+i)&3));
		}
		out.write('\n');
	}

	private static void writeFourDigits(BufferedWriter out, int value) throws IOException{
		out.write((char)('0'+((value/1000)%10)));
		out.write((char)('0'+((value/100)%10)));
		out.write((char)('0'+((value/10)%10)));
		out.write((char)('0'+(value%10)));
	}

	private static void checkFasta(String label, String workload, File fasta, int expected,
			int expectedLength) throws IOException{
		final String expectedSecond=workload.equals("mixed") ? "eligible_" : "all_";
		int headers=0, sequences=0, badLength=0, badBase=0, badSecond=0, seqLength=0;
		try(BufferedReader in=new BufferedReader(new InputStreamReader(
				new FileInputStream(fasta), StandardCharsets.US_ASCII))){
			for(String line; (line=in.readLine())!=null; ){
				if(line.startsWith(">")){
					headers++;
					if(!hasExpectedSecondParent(line, expectedSecond)){badSecond++;}
					if(seqLength>0){
						sequences++;
						if(seqLength!=expectedLength){badLength++;}
						seqLength=0;
					}
				}else{
					seqLength+=line.length();
					for(int i=0; i<line.length(); i++){
						final char c=line.charAt(i);
						if(c!='A' && c!='C' && c!='G' && c!='T'){badBase++;}
					}
				}
			}
		}
		if(seqLength>0){
			sequences++;
			if(seqLength!=expectedLength){badLength++;}
		}
		System.out.println(label + '\t' + workload + '\t' + headers + '\t' + sequences
				+ '\t' + badLength + '\t' + badBase + '\t' + badSecond);
		if(headers!=expected || sequences!=expected || badLength!=0 || badBase!=0
				|| badSecond!=0){
			throw new IllegalStateException("Invalid FASTA for " + label + '/' + workload);
		}
	}

	private static boolean hasExpectedSecondParent(String line, String expectedSecond){
		final int second=line.lastIndexOf(" ~ ");
		if(second<0){return false;}
		final int start=second+3;
		if(!line.startsWith(expectedSecond, start)){return false;}
		if(line.length()!=start+expectedSecond.length()+4){return false;}
		for(int i=start+expectedSecond.length(); i<line.length(); i++){
			final char c=line.charAt(i);
			if(c<'0' || c>'9'){return false;}
		}
		return true;
	}
}
