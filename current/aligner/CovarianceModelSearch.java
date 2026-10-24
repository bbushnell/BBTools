package aligner;

import java.util.Locale;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.ObjectSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Exact best root interval per strand. This is raw global-CM CYK, not a
 * biological caller: no score threshold, null3, calibrated E-value, local model
 * ends, truncation, bands, filters or multi-hit selection. Forward target bounds
 * are one-based inclusive; parsetrees retain oriented input coordinates.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelSearch {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CovarianceModelSearch.class, false);final Parser p=new Parser();
		String cm=null, trace=null, coordinates=null, ties=null;long maxCells=200000000;int maxWindow=0;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("model")){cm=value;}else if(key.equals("trace")){trace=value;}else if(key.equals("coordinates")){coordinates=value;}
			else if(key.equals("ties")){ties=value;}else if(key.equals("maxcells")){maxCells=Parse.parseKMG(value);}
			else if(key.equals("maxwindow")){maxWindow=Parse.parseIntKMG(value);}
			else{require(p.parse(arg, key, value), "Unknown argument: "+key);}}
		require(cm!=null && p.in1!=null && p.out1!=null && trace!=null && coordinates!=null && ties!=null && maxWindow>0,
			"Bind model, in, out, trace, coordinates, ties and positive maxwindow; oracle W must match");
		final String[] paths={cm,p.in1,p.out1,trace,coordinates,ties};for(int i=0; i<paths.length; i++){for(int j=0; j<i; j++){require(!paths[i].equals(paths[j]), "Input/output paths must differ");}}
		CovarianceModelAlphabet.configureStandaloneInput();
		Shared.setThreads(1);final CovarianceModel model=CovarianceModelParser.read(cm);final CovarianceModelCyk scorer=new CovarianceModelCyk(model, maxCells);
		final ByteStreamWriter scores=writer(p.out1, p.overwrite), traces=writer(trace, p.overwrite), coords=writer(coordinates, p.overwrite), tied=writer(ties, p.overwrite);
		final Streamer reader=StreamerFactory.makeStreamer(FileFormat.testInput(p.in1, FileFormat.FASTA, null, true, true), null, true, -1, false, true, 1);
		final ByteBuilder b=new ByteBuilder();final ObjectSet<String> ids=new ObjectSet<String>(String.class);int count=0, unsupported=0, impossible=0;boolean success=false;
		try{scores.println("sequence\tstrand\tmodel\tlength\tmaxWindow\tstatus\trawCykBits\ttargetStart1\ttargetStop1\tbestIntervalTies\tmatchedResidues\tinsertedResidues\tdeletedModelPositions\tmodelFrom1\tmodelTo1\tcells\tevalue\tevalueReason");
			traces.println("sequence\tstrand\tstep\tstate\ttype\torientedFrom1\torientedTo1\tchoice\tleftChildStep\trightChildStep");
			coords.println("sequence\tstrand\ttargetPosition1\torientedPosition1\tstate\ttype\tmodelPosition1\tinsertAfterModelPosition");
			tied.println("sequence\tstrand\ttieIndex\ttargetStart1\ttargetStop1\trawCykBits");reader.start();
			for(ListNum<Read> batch; (batch=reader.nextList())!=null;){for(Read r:batch.list){
				final String id=first(r.id);require(ids.add(id), "Duplicate FASTA ID makes strand/coordinate joins ambiguous: "+id);count++;
				for(int orientation=0; orientation<2; orientation++){final boolean reverse=orientation==1;final char strand=reverse ? '-' : '+';final CovarianceModelCyk.Result result;
					try{result=scorer.search(reverse ? reverseComplement(r.bases) : r.bases, maxWindow);}
					catch(CovarianceModelCyk.UnsupportedSequenceException e){unsupported++;
						scores.print(b.clear().append(id).tab().append(strand).tab().append(model.name).tab().append(r.length()).tab().append(maxWindow)
							.append("\tUNSUPPORTED_SEQUENCE\tNA\t0\t0\t0\tNA\tNA\tNA\tNA\tNA\t0\tNA\tuncalibrated_mode\n"));
						pp.outstream.println(id+" "+strand+": "+e.getMessage());continue;}
					final boolean finite=Float.isFinite(result.score);if(!finite){impossible++;}
					final int start=finite ? forwardPosition(reverse ? result.targetTo : result.targetFrom, r.length(), reverse) : 0;
					final int stop=finite ? forwardPosition(reverse ? result.targetFrom : result.targetTo, r.length(), reverse) : 0;
					scores.print(b.clear().append(id).tab().append(strand).tab().append(model.name).tab().append(r.length()).tab().append(maxWindow).tab().append(finite ? "OK" : "IMPOSSIBLE").tab().appendSlow(result.score)
						.tab().append(start).tab().append(stop).tab().append(result.rootStarts.size).tab().append(result.matched).tab().append(result.inserted).tab().append(result.deleted)
						.tab().append(result.modelFrom).tab().append(result.modelTo).tab().append(result.cells).append("\tNA\tuncalibrated_mode\n"));
					for(int n=0; n<result.rootStarts.size; n++){
						final int a=forwardPosition(reverse ? result.rootStops.get(n) : result.rootStarts.get(n), r.length(), reverse), z=forwardPosition(reverse ? result.rootStarts.get(n) : result.rootStops.get(n), r.length(), reverse);
						tied.print(b.clear().append(id).tab().append(strand).tab().append(n).tab().append(a).tab().append(z).tab().appendSlow(result.score).nl());}
					for(int n=0; n<result.states.size; n++){final int v=result.states.get(n);traces.print(b.clear().append(id).tab().append(strand).tab().append(n).tab().append(v).tab().append(CovarianceModel.STATE_NAMES[model.type[v]])
						.tab().append(result.from.get(n)).tab().append(result.to.get(n)).tab().append(result.choice.get(n)).tab().append(result.left.get(n)).tab().append(result.right.get(n)).nl());}
					if(finite){for(int i=result.targetFrom; i<=result.targetTo; i++){final int v=result.stateAtBase[i-1];coords.print(b.clear().append(id).tab().append(strand).tab().append(forwardPosition(i, r.length(), reverse)).tab().append(i)
						.tab().append(v).tab().append(CovarianceModel.STATE_NAMES[model.type[v]]).tab().append(result.modelPosition[i-1]).tab().append(result.insertAnchor[i-1]).nl());}}
				}
			}reader.returnList(batch);}
			require(!reader.errorState(), "FASTA streaming failed; search output is incomplete");success=true;
		}finally{if(!success){reader.close();}final boolean a=scores.poisonAndWait(), z=traces.poisonAndWait(), c=coords.poisonAndWait(), t=tied.poisonAndWait();require(!(a || z || c || t), "Search output failed");}
		pp.outstream.println("CM_SEARCH_PASS sequences="+count+" strands="+(2L*count)+" unsupported="+unsupported+" impossible="+impossible);Shared.closeStream(pp.outstream);
	}

	public static int forwardPosition(int oriented, int length, boolean reverse){
		require(oriented>=1 && oriented<=length, "Emitted coordinate must be inside the input before strand conversion");return reverse ? length-oriented+1 : oriented;
	}
	public static byte[] reverseComplement(byte[] input){
		assert(input!=null):"Strand conversion cannot invent missing input";
		final byte[] result=new byte[input.length];for(int i=0; i<input.length; i++){
			result[input.length-i-1]=CovarianceModelAlphabet.complement(input[i]);}return result;
	}
	private static String first(String s){require(s!=null && !s.isEmpty(), "FASTA IDs must be nonempty");int n=0;while(n<s.length() && !Character.isWhitespace(s.charAt(n))){n++;}return s.substring(0, n);}
	private static ByteStreamWriter writer(String path, boolean overwrite){final ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(path, FileFormat.TEXT, null, true, overwrite, false, false));w.start();return w;}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
}
