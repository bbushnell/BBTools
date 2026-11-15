package test.icecream;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;

import fileIO.FileFormat;
import icecream.ZMW;
import icecream.ZMWStreamer;
import shared.Shared;
import stream.Read;

/** Checks complete ZMW delivery against literal well IDs and read counts.
 * @author Brian Bushnell, Sucrose
 */
public final class ZMWStreamerLimitTest{

	/** Runs finite EOF/limit cases through the production FASTA or SAM input path.
	 * @param args Required outdir=new-directory; optional backend=cris/streamer/both and scope=limits/movies/all
	 * @throws Exception If fixture input creation or producer joining fails
	 */
	public static void main(String[] args) throws Exception{
		Path outdir=null;
		String backend="both";
		String scope="limits";
		for(String arg : args){
			if(arg.startsWith("outdir=")){outdir=Paths.get(arg.substring(7));}
			else if(arg.startsWith("backend=")){backend=arg.substring(8);}
			else if(arg.startsWith("scope=")){scope=arg.substring(6);}
			else{throw new IllegalArgumentException("Unknown fixture argument: "+arg);}
		}
		require(outdir!=null, "Specify outdir=<new-directory> to preserve each finite fixture input");
		require(backend.equals("cris") || backend.equals("streamer") || backend.equals("both"),
				"backend must be cris, streamer, or both: "+backend);
		require(scope.equals("limits") || scope.equals("movies") || scope.equals("all"),
				"scope must be limits, movies, or all: "+scope);
		require(ZMWStreamer.class.desiredAssertionStatus(), "Production assertions must be enabled");
		Files.createDirectory(outdir);
		Shared.setThreads(1);
		final boolean originalStreamer=ZMWStreamer.useStreamer;
		int failures=0;
		try{
			System.out.println("backend\tcase\tlimit\texpected\tobserved\tresult");
			if(!scope.equals("movies")){
				if(!backend.equals("streamer")){failures+=runCases(outdir, false);}
				if(!backend.equals("cris")){failures+=runCases(outdir, true);}
			}
			if(!scope.equals("limits")){
				if(!backend.equals("streamer")){failures+=runMovieCases(outdir, false);}
				if(!backend.equals("cris")){failures+=runMovieCases(outdir, true);}
			}
		}finally{
			ZMWStreamer.useStreamer=originalStreamer;
		}
		require(failures==0, failures+" ZMW limit cases differed from the literal well/count oracle; see all case rows");
	}

	/** Keeps expectations independent of the production buffering and limit conditions. */
	private static int runCases(Path outdir, boolean sam) throws Exception{
		int failures=0;
		failures+=observe(outdir, sam, "eof_below_limit", 2, 2, "[11:2]");
		failures+=observe(outdir, sam, "eof_at_limit", 3, 2, "[11:2,22:1]");
		failures+=observe(outdir, sam, "input_above_limit", 5, 2, "[11:2,22:1]");
		failures+=observe(outdir, sam, "unlimited", 5, -1, "[11:2,22:1,33:2]");
		return failures;
	}

	/** Consumes the real queue through its terminal, then checks well contents and producer count. */
	private static int observe(Path outdir, boolean sam, String label, int readCount,
			long limit, String expected) throws Exception{
		final String backend=sam ? "streamer" : "cris";
		final Path input=outdir.resolve(backend+"-"+label+(sam ? ".sam" : ".fa"));
		writeInput(input, sam, readCount);
		ZMWStreamer.useStreamer=sam;
		final FileFormat ff=FileFormat.testInput(input.toString(), sam ? FileFormat.SAM : FileFormat.FASTA,
				null, true, true);
		final ZMWStreamer streamer=new ZMWStreamer(ff, 4, -1, limit);
		final Thread producer=streamer.runStreamer(true);
		final ArrayList<ZMW> wells=new ArrayList<ZMW>();
		for(ZMW well=streamer.nextZMW(); well!=null; well=streamer.nextZMW()){
			wells.add(well);
		}
		producer.join();
		require(!streamer.errorState, backend+" input failed in finite case "+label);
		final String observed=describe(wells);
		final boolean pass=expected.equals(observed) && streamer.ZMWs==wells.size();
		System.out.println(backend+"\t"+label+"\t"+limit+"\t"+expected+"\t"+observed+"\t"+(pass ? "PASS" : "FAIL"));
		return pass ? 0 : 1;
	}

	/** Distinct movies must remain separate even when their numeric well IDs agree. */
	private static int runMovieCases(Path outdir, boolean sam) throws Exception{
		return observeMovie(outdir, sam, "same_number_same_movie",
				new String[]{"movieA/42/0_4", "movieA/42/4_8"}, "[movieA/42:2]")+
				observeMovie(outdir, sam, "same_number_different_movie",
						new String[]{"movieA/42/0_4", "movieB/42/0_4"}, "[movieA/42:1,movieB/42:1]")+
				observeMovie(outdir, sam, "unparsed_numeric_fallback",
						new String[]{"fallback/42/noCoords", "fallback/42/noCoords"}, "[fallback/42:1,fallback/42:1]");
	}

	/** Two literal records exercise movie identity without rerunning the accepted output-limit cases. */
	private static int observeMovie(Path outdir, boolean sam, String label, String[] ids, String expected) throws Exception{
		final String backend=sam ? "streamer" : "cris";
		final String[] bases={"ACGT", "TGCA"};
		final Path input=outdir.resolve(backend+"-"+label+(sam ? ".sam" : ".fa"));
		final StringBuilder data=new StringBuilder();
		if(sam){data.append("@HD\tVN:1.6\tSO:unknown\n");}
		for(int i=0; i<2; i++){
			if(sam){data.append(ids[i]).append("\t4\t*\t0\t0\t*\t*\t0\t0\t").append(bases[i]).append("\tIIII\n");}
			else{data.append('>').append(ids[i]).append('\n').append(bases[i]).append('\n');}
		}
		Files.write(input, data.toString().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		ZMWStreamer.useStreamer=sam;
		final FileFormat ff=FileFormat.testInput(input.toString(), sam ? FileFormat.SAM : FileFormat.FASTA, null, true, true);
		final ZMWStreamer streamer=new ZMWStreamer(ff, 4, -1, -1);
		final Thread producer=streamer.runStreamer(true);
		final StringBuilder groups=new StringBuilder("[");
		int reads=0, wells=0;
		boolean contents=true;
		for(ZMW well=streamer.nextZMW(); well!=null; well=streamer.nextZMW()){
			require(!well.isEmpty(), "A movie group must contain at least one literal read");
			if(wells++>0){groups.append(',');}
			final String first=well.get(0).id;
			groups.append(first, 0, first.lastIndexOf('/')).append(':').append(well.size());
			for(Read read : well){
				contents&=reads<2 && ids[reads].equals(read.id) &&
						bases[reads].equals(new String(read.bases, StandardCharsets.US_ASCII));
				reads++;
			}
		}
		producer.join();
		require(!streamer.errorState, backend+" input failed in movie case "+label);
		final String observed=groups.append(']').toString();
		final boolean pass=contents && reads==2 && streamer.ZMWs==wells && expected.equals(observed);
		System.out.println(backend+"\t"+label+"\t-1\t"+expected+"\t"+observed+"\t"+(pass ? "PASS" : "FAIL"));
		return pass ? 0 : 1;
	}

	/** Writes at most five literal reads; the two input formats have identical names and bases. */
	private static void writeInput(Path path, boolean sam, int readCount) throws Exception{
		require(readCount==2 || readCount==3 || readCount==5, "Finite input must end at a complete well boundary: "+readCount);
		final StringBuilder text=new StringBuilder();
		if(sam){text.append("@HD\tVN:1.6\tSO:unknown\n");}
		for(int i=0; i<readCount; i++){
			if(sam){text.append(IDS[i]).append("\t4\t*\t0\t0\t*\t*\t0\t0\t").append(BASES[i]).append("\tIIII\n");}
			else{text.append('>').append(IDS[i]).append('\n').append(BASES[i]).append('\n');}
		}
		Files.write(path, text.toString().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Checks literal read order/bases as well as reporting the well/count oracle independently of PBHeader. */
	private static String describe(ArrayList<ZMW> wells){
		final StringBuilder text=new StringBuilder("[");
		int readIndex=0;
		for(int i=0; i<wells.size(); i++){
			final ZMW well=wells.get(i);
			require(!well.isEmpty(), "Production must not deliver an empty nonterminal well");
			if(i>0){text.append(',');}
			final String first=well.get(0).id;
			final int slash=first.indexOf('/', "movie/".length());
			require(first.startsWith("movie/") && slash>6, "Unexpected fixture header: "+first);
			final String prefix=first.substring(0, slash+1);
			text.append(first, 6, slash).append(':').append(well.size());
			for(Read read : well){
				require(readIndex<IDS.length && read.id.equals(IDS[readIndex]) && read.id.startsWith(prefix),
						"Well delivery must preserve complete, ordered literal read IDs: "+read.id+", index="+readIndex);
				require(BASES[readIndex].equals(new String(read.bases, StandardCharsets.US_ASCII)),
						"Well grouping must preserve sequence bytes for "+read.id);
				readIndex++;
			}
		}
		return text.append(']').toString();
	}

	/** Keeps oracle and setup failures visible independently of Java assertion settings. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Literal finite input; expected groups above are specified separately. */
	private static final String[] IDS={"movie/11/0_4", "movie/11/4_8", "movie/22/0_4", "movie/33/0_4", "movie/33/4_8"};
	private static final String[] BASES={"ACGT", "TGCA", "AAAA", "CCCC", "GGGG"};

	/** Static-only fixture. */
	private ZMWStreamerLimitTest(){}
}
