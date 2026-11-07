package test.prok;

import java.io.File;

import parse.Parse;
import prok.NcrnaOccurrenceCheck;
import prok.NcrnaTextIO.Input;
import structures.ByteBuilder;

/** Checks actual sequence parity and rejects independently corrupted finite inputs.
 * Reuses the literal native-extraction fixture, not a biological corpus.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaOccurrenceCheckTest{

	public static void main(final String[] args){
		String out=null;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0], b=split.length>1 ? split[1] : null;
			if(a.equalsIgnoreCase("out")){out=b;}
			else{throw new IllegalArgumentException("Unknown fixture parameter: "+arg);}
		}
		if(out==null || !new File(out).mkdir()){throw new IllegalArgumentException("Fresh fixture directory required");}
		final String base=new File(out, "base").getAbsolutePath();
		NcrnaOccurrenceInputsTest.main(new String[]{"out="+base});
		final String ledger=base+"/prepared/occurrence_loci.tsv", contexts=base+"/contexts.fa", cores=base+"/cores.fa";
		final String summary=new File(out, "check.tsv").getAbsolutePath();
		invoke(ledger, contexts, cores, summary);
		if(!text(summary).endsWith("2\t2\t2\t51\t41\t1\n")){throw new AssertionError("Literal complete sequence/base/contig counts");}
		final String bad=new File(out, "bad.fa").getAbsolutePath(), failure=new File(out, "failure.tsv").getAbsolutePath();
		final String coreText=text(cores), contextText=text(contexts);
		write(bad, coreText.replaceFirst("AAAAAAAAAAAAAAAAAAAA", "TAAAAAAAAAAAAAAAAAAA"));
		fails(ledger, contexts, bad, failure, "Core differs");
		write(bad, coreText.replaceFirst("AAAAAAAAAAAAAAAAAAAA", "aaaaaaaaaaaaaaaaaaaa"));
		fails(ledger, contexts, bad, failure, "Core differs");
		write(bad, coreText.replace("ID=occ2", "ID=occ1"));
		fails(ledger, contexts, bad, failure, "Duplicate core");
		write(bad, contextText.replace("lflank=0", "lflank=1"));
		fails(ledger, bad, cores, failure, "left flank");
		write(bad, contextText.replace("source=gA", "source=gB"));
		fails(ledger, bad, cores, failure, "Source identity");
		write(bad, contextText.replaceAll("contig_gc=[0-9.]+", "contig_gc=NaN"));
		fails(ledger, bad, cores, failure, "GC must be finite");
		write(bad, contextText.replaceFirst("contig_gc=[0-9.]+", "contig_gc=0.5"));
		fails(ledger, bad, cores, failure, "inconsistent caller GC");
		write(bad, contextText.replace("ID=occ2", "ID=unknown"));
		fails(ledger, bad, cores, failure, "Sequence ID is absent");
		write(bad, contextText+contextText);
		fails(ledger, bad, cores, failure, "Context FASTA count");
		System.err.println("NCRNA_OCCURRENCE_CHECK_FINITE_PASS");
	}

	private static void invoke(final String ledger, final String contexts, final String cores, final String out){
		NcrnaOccurrenceCheck.main(new String[]{"in="+ledger, "contexts="+contexts, "cores="+cores, "out="+out, "expected=2"});
	}
	private static void fails(final String ledger, final String contexts, final String cores, final String out, final String reason){
		boolean failed=false;
		try{invoke(ledger, contexts, cores, out);}
		catch(IllegalArgumentException expected){
			if(!expected.getMessage().contains(reason)){throw new AssertionError("Unexpected failure: "+expected);}
			failed=true;
		}
		if(!failed || new File(out).exists()){throw new AssertionError("Corrupt inputs must fail before output: "+reason);}
	}
	private static String text(final String path){
		final ByteBuilder buffer=new ByteBuilder();
		try(Input input=new Input(path)){
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){buffer.append(row).nl();}
		}
		return buffer.toString();
	}
	private static void write(final String path, final String text){NcrnaLearningPopulationTest.write(path, text);}
}
