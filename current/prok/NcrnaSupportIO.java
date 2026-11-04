package prok;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import prok.NcrnaSupport.Columns;
import prok.NcrnaSupport.Graph;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Native byte-oriented reader/writer for the shared scored support graph.
 * Version 2 carries explicit ORIGINAL/REDUCED scope and dense named model axes.
 * A48 stores float bits without threshold-changing decimal rounding. Reads are
 * validated in one pass; writing never hashes or reopens its output.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSupportIO{

	/*--------------------------------------------------------------*/
	/*----------------            Reading           ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads one complete graph, failing on malformed framing or inconsistent edges. */
	public static Graph read(final String path){
		try(Input owner=new Input(path)){
			final GraphReader input=new GraphReader(owner.reader);
			final LineParser1 lp=input.parser;
			input.header("#ncrna_support_graph", 2);
			if(integerField(lp, 1)!=2){throw bad("Unsupported support graph version");}
			input.header("#scope", 2);
			final boolean reduced=lp.termEquals("REDUCED", 1);
			if(!reduced && !lp.termEquals("ORIGINAL", 1)){throw bad("Graph scope must be ORIGINAL or REDUCED");}
			input.header("#aligner", 2);
			if(!lp.termEquals("quantum", 1)){throw bad("This graph contract requires Quantum identities");}
			input.header("#cutoff_a48", 2);
			final float cutoff=floatField(lp, 1);
			input.header("#models", 2);
			final int modelCount=positiveCount(lp, 1);
			final String[] models=new String[modelCount];
			for(int model=0; model<modelCount; model++){
				input.header("#model", 3);
				if(integerField(lp, 1)!=model){throw bad("Model axis must contain each dense index in order");}
				models[model]=lp.parseString(2);
			}
			input.header("#queries", 2);
			final int queries=positiveCount(lp, 1);
			input.header("query_id", 4);
			if(!lp.termEquals("best_model_index", 1) || !lp.termEquals("best_identity_a48", 2) || !lp.termEquals("supporters", 3)){
				throw bad("Unexpected support graph data columns");
			}
			final Columns data=new Columns();
			for(int query=0; query<queries; query++){
				input.row(4);
				final String name=lp.parseString(0);
				final int best=lp.termEquals("NONE", 1) ? -1 : integerField(lp, 1);
				final float identity=lp.termEquals("NA", 2) ? Float.NaN : floatField(lp, 2);
				if(!lp.termEquals(".", 3)){
					lp.setBounds(3);
					parseEdges(lp.line(), lp.a(), lp.b(), data);
				}
				data.finish(name, best, identity);
			}
			if(input.file.nextLine()!=null){throw bad("Trailing data after the declared support queries");}
			return data.graph(models, cutoff, reduced);
		}
	}

	/** Parses model:identity tokens using only the first colon; A48 may contain colons. */
	private static void parseEdges(final byte[] line, final int from, final int to, final Columns data){
		if(from==to){throw bad("Empty supporter field; use a dot for no support");}
		int begin=from;
		for(int end=from; end<=to; end++){
			if(end<to && line[end]!=','){continue;}
			int colon=begin;
			while(colon<end && line[colon]!=':'){colon++;}
			if(colon==begin || colon==end){throw bad("Support entries must be modelIndex:A48Identity");}
			data.models.add(NcrnaNumbers.unsignedInt(line, begin, colon));
			data.identities.add(decodeFloat(line, colon+1, end));
			begin=end+1;
		}
	}

	private static int positiveCount(final LineParser1 lp, final int term){
		final int count=integerField(lp, term);
		if(count<1 || count>=Integer.MAX_VALUE){throw bad("Graph axes require positive, representable counts");}
		return count;
	}

	private static int integerField(final LineParser1 lp, final int term){
		lp.setBounds(term);
		return NcrnaNumbers.unsignedInt(lp.line(), lp.a(), lp.b());
	}

	private static float floatField(final LineParser1 lp, final int term){
		lp.setBounds(term);
		return decodeFloat(lp.line(), lp.a(), lp.b());
	}

	/** Validates the native unsigned float word before calling the existing decoder. */
	private static float decodeFloat(final byte[] line, final int from, final int to){
		if(to<=from || to-from>6){throw bad("A48 float requires one to six digits");}
		for(int i=from; i<to; i++){
			if(line[i]<'0' || line[i]>'o'){throw bad("Invalid A48 float digit");}
		}
		final long bits=Parse.parseLongA48(line, from, to);
		if(bits>0xffffffffL){throw bad("A48 word exceeds the 32-bit float representation");}
		final float value=Float.intBitsToFloat((int)bits);
		if(!(value>=0 && value<=1)){throw bad("Support identity must be finite and in [0,1]");}
		return value;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Writing           ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes a new graph through ByteStreamWriter; compression follows its suffix. */
	public static void write(final Graph graph, final String path){
		assert(graph!=null) : "Support serialization needs the graph's named axes and scope";
		try(Output owner=new Output(path)){
			final GraphWriter output=new GraphWriter(owner.writer);
			final ByteBuilder bb=output.buffer;
			bb.append("#ncrna_support_graph\t2\n#scope\t").append(graph.reduced ? "REDUCED" : "ORIGINAL").nl();
			bb.append("#aligner\tquantum\n#cutoff_a48\t").appendFloatA48(graph.cutoff).nl();
			bb.append("#models\t").append(graph.modelNames.length).nl();
			output.flush();
			for(int model=0; model<graph.modelNames.length; model++){
				bb.append("#model\t").append(model).tab().append(graph.modelNames[model]).nl();
				output.flush();
			}
			bb.append("#queries\t").append(graph.queryNames.length).nl();
			bb.append("query_id\tbest_model_index\tbest_identity_a48\tsupporters\n");
			output.flush();
			for(int query=0; query<graph.queryNames.length; query++){
				bb.append(graph.queryNames[query]).tab();
				if(graph.bestModels[query]<0){
					bb.append("NONE\tNA");
				}else{
					bb.append(graph.bestModels[query]).tab().appendFloatA48(graph.bestIdentities[query]);
				}
				bb.tab();
				final int from=graph.starts[query], to=graph.starts[query+1];
				if(from==to){bb.append('.');}
				for(int edge=from; edge<to; edge++){
					if(edge>from){bb.append(',');}
					bb.append(graph.models[edge]).append(':').appendFloatA48(graph.identities[edge]);
				}
				bb.nl();
				output.flush();
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Format Helpers       ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrows an already-owned reader and retains reusable graph field offsets. */
	private static final class GraphReader{
		GraphReader(final ByteFile file_){file=file_;}

		void header(final String name, final int terms){
			row(terms);
			if(!parser.termEquals(name, 0)){throw bad("Expected graph header "+name);}
		}

		void row(final int terms){
			final byte[] line=file.nextLine();
			if(line==null){throw bad("Truncated support graph");}
			parser.set(line);
			if(parser.terms()!=terms){throw bad("Expected "+terms+" tab-separated graph fields");}
		}

		final ByteFile file;
		final LineParser1 parser=new LineParser1('\t');
	}

	/** Borrows an already-owned writer and retains one reusable record buffer. */
	private static final class GraphWriter{
		GraphWriter(final ByteStreamWriter writer_){writer=writer_;}

		void flush(){
			writer.print(buffer);
			buffer.clear();
		}

		final ByteStreamWriter writer;
		final ByteBuilder buffer=new ByteBuilder();
	}

	private static IllegalArgumentException bad(final String message){return new IllegalArgumentException(message);}
	private NcrnaSupportIO(){}
}
