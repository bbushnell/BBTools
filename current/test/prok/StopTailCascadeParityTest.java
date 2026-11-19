package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.StopTailCascadeResources;
import prot.DigestSuffix;
import prot.FamilyScoringPackage;
import prot.FamilyTailCascade;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Writes every assignment field for a prebound real-tail panel through either
 * evaluator or caller loading. Run both routes separately and compare exact bytes.
 * @author Keqing
 */
public final class StopTailCascadeParityTest{

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, StopTailCascadeParityTest.class, false);
		try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(String[] args) throws Exception{
		final Parser parser=new Parser();
		String route=null, binding=null, bindingPin=null, inputPin=null;
		int families=-1, queries=-1;
		double scale=-1;
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
			if(key.equals("route")){route=value;}
			else if(key.equals("binding")){binding=value;}
			else if(key.equals("bindingsha80")){bindingPin=value;}
			else if(key.equals("insha80")){inputPin=value;}
			else if(key.equals("families")){families=Parse.parseIntKMG(value);}
			else if(key.equals("queries")){queries=Parse.parseIntKMG(value);}
			else if(key.equals("marginscale")){scale=Double.parseDouble(value);}
			else if(!parser.parse(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		require("evaluator".equals(route) || "caller".equals(route), "route must be evaluator or caller");
		require(binding!=null && parser.in1!=null && parser.out1!=null, "Require binding, in and out");
		require(parser.in2==null && parser.out2==null && !parser.append && parser.maxReads<0,
			"Parity requires the entire prebound unpaired panel, without append or read limits");
		require(families>0 && queries>0 && Double.isFinite(scale) && scale>=0,
			"Require positive expected families/queries and a finite nonnegative margin scale");
		final Path input=Paths.get(parser.in1).toRealPath(), output=Paths.get(parser.out1).toAbsolutePath();
		require(!Files.exists(output), "Parity output must be new: "+output);
		verify(input.toString(), inputPin); verify(binding, bindingPin);
		Shared.AMINO_IN=true;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final int threads=parser.workers()>0 ? parser.workers() : 1;
		final FamilyScoringPackage evaluator;
		final FamilyTailCascade.LoadedClassifier classifier;
		if(route.equals("evaluator")){
			evaluator=new FamilyScoringPackage(Paths.get(binding));
			classifier=evaluator.load(threads).withEvaluationMarginScale(scale);
		}else{
			evaluator=null;
			final StopTailCascadeResources caller=StopTailCascadeResources.read(binding);
			require(Double.compare(scale, caller.marginScale)==0, "Caller manifest must bind the preselected margin exactly");
			caller.protectOutputs(new String[]{output.toString()});
			classifier=caller.load(threads);
		}
		require(classifier.familyCount()==families && !classifier.usesEarlyRejection(),
			"Parity must use the expected live family count and explicit exhaustive index policy");
		final FamilyTailCascade.Worker worker=classifier.newWorker();
		final Streamer stream=StreamerFactory.makeStreamer(FileFormat.testInput(input.toString(), FileFormat.FASTA,
			null, true, true), null, true, -1, false, true, 1);
		final ByteStreamWriter out=new ByteStreamWriter(output.toString(), false, false, false);
		final ByteBuilder row=new ByteBuilder();
		long count=0, accepted=0, rejected=0;
		String firstAccepted=null, firstRejected=null;
		stream.start(); out.start();
		try{
			out.print("ordinal\tquery\tassigned\tfamily\tterminal\troster_rank\tddl_matches\tddl_position\t"
				+"dp_attempt\tdp_attempts\tscore\tidentity\tblosum\tscore_measured\tidentity_measured\t"
				+"blosum_measured\n");
			for(ListNum<Read> batch=stream.nextList();batch!=null;batch=stream.nextList()){
				for(Read read:batch.list){
					require(read!=null && read.mate==null && read.bases!=null && read.bases.length==200,
						"Parity panel requires unpaired exactly200-aa real tails");
					final String id=Tools.trimToWhitespace(read.id);
					require(id!=null && !id.isEmpty(), "Real-tail parity requires a nonempty FASTA identifier");
					final FamilyTailCascade.Assignment a=worker.assign(id, read.bases);
					row.clear().append(count++).tab().append(id).tab().append(a.assigned()).tab()
						.append(a.familyId==null ? "NA" : a.familyId).tab().append(a.terminal).tab()
						.append(a.rosterRank).tab().append(a.ddlMatches).tab().append(a.ddlPosition).tab()
						.append(a.dpAttempt).tab().append(a.dpAttempts).tab().append(Double.toString(a.score)).tab()
						.append(Float.toString(a.identity)).tab().append(Float.toString(a.blosum)).tab()
						.append(a.scoreMeasured).tab().append(a.identityMeasured).tab().append(a.blosumMeasured).nl();
					out.print(row);
					if(a.assigned()){accepted++; if(firstAccepted==null){firstAccepted=id;}}
					else{rejected++; if(firstRejected==null){firstRejected=id;}}
				}
			}
		}finally{
			final boolean inputFailed=ReadWrite.closeStream(stream), outputFailed=out.poisonAndWait();
			require(!inputFailed && !outputFailed, "Parity input or assignment output failed");
		}
		require(count==queries && accepted>0 && rejected>0,
			"Parity must cover all expected queries and both outcomes: "+count+"/"+accepted+"/"+rejected);
		if(evaluator!=null){evaluator.verify();}
		verify(binding, bindingPin); verify(input.toString(), inputPin);
		System.out.println("STOP_TAIL_PARITY_PASS route="+route+" queries="+count
			+" accepted="+accepted+" rejected="+rejected
			+" first_accepted="+firstAccepted+" first_rejected="+firstRejected);
	}

	private static void verify(String path, String pin){
		require(pin!=null && pin.matches("[0-9a-f]{20}") && pin.equals(DigestSuffix.file(path)),
			"Input binding differs: "+path);
	}
	private static void require(boolean valid, String why){if(!valid){throw new IllegalArgumentException(why);}}
}
