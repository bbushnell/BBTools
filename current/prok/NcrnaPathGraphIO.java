package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.ByteFile;
import map.LongHashSet;
import parse.LineParser1;
import structures.ByteBuilder;
import structures.FloatList;
import structures.IntList;
import structures.LongList;

/** Versioned, contig-framed path metadata; DNA never enters this format.
 * The writer reuses one contig buffer. Only contigs with a returned target RNA
 * retain candidates; all other contigs retain their exact baseline report rows.
 * @author Brian, Raiden
 */
public final class NcrnaPathGraphIO{
	public interface Policies{
		int index(String source);
		GeneCallerPathPolicy policy(int source);
		default void verifyScoreBasis(String family, String model, int nativeModel){
			throw new IllegalStateException("Bound graph reading requires actual family/model axis validation");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Native capture         ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Capture implements GeneCallerPathTrace.Observer{
		public Capture(final String source_, final GeneCaller caller_){
			require(source_!=null && !source_.isEmpty() && caller_!=null, "Capture requires a declared source and its native caller");
			source=source_;
			caller=caller_;
		}
		public void start(final long read_, final String contig_, final int length_){
			require(!active && read_>=0 && contig_!=null && !contig_.isEmpty() && length_>=0,
				"Start one source-qualified physical contig after consuming the previous record");
			read=read_;
			contig=contig_;
			length=length_;
			active=true;
			begun=ended=false;
			targets=basisRows=0;
			body.clear();
			output.clear();
		}
		@Override
		public void begin(final String name, final int contigLength, final int candidates){
			require(active && !begun && contigLength==length && candidates>=0
				&& (candidates==0 ? name==null : contig.equals(name)), "Native graph must match the current read identity and length");
			if(policy==null){policy=caller.pathPolicy();}else{policy.validate(caller);}
			count=candidates;
			begun=true;
			policy.append(body).nl();
		}
		@Override
		public boolean candidate(long id, int start, int stop, int strand, int frame, int type,
				int stopLeft, float score, float bonus, boolean nativeMeasured, float nativeScore){
			body.append("C\t").append(id).tab().append(start).tab().append(stop).tab().append(strand)
				.tab().append(frame).tab().append(type).tab().append(stopLeft).tab().append(bits(score))
				.tab().append(bits(bonus)).tab().append(nativeMeasured ? 1 : 0).tab().append(bits(nativeScore)).nl();
			return false;
		}
		@Override
		public void nodeClass(final long id, final String name){
			final int kind=Orf.class.getName().equals(name) ? 0 : NcrnaReportedOrf.class.getName().equals(name) ? 1 : -1;
			require(kind>=0, "Saved paths do not support unrecorded candidate-subtype state: "+name);
			body.append("K\t").append(id).tab().append(kind).nl();
		}
		@Override
		public void metadata(long id, int flipped, int startCodon, int stopCodon,
				boolean startTruncated, boolean stopTruncated, float startScore, float stopScore,
				float kmerScore, float oppositePenalty, String model, String anticodon, String family,
				int partial, float acceptedRrnaIdentity, boolean deferred, int reportStart, int reportStop,
				boolean proteinEvidence, boolean ssuEvidence){
			require(!proteinEvidence && !ssuEvidence && oppositePenalty==0, "Capture excludes unsupported evidence-dependent policies");
			body.append("M\t").append(id).tab().append(flipped).tab().append(startCodon).tab().append(stopCodon)
				.tab().append(startTruncated ? 1 : 0).tab().append(stopTruncated ? 1 : 0)
				.tab().append(bits(startScore)).tab().append(bits(stopScore)).tab().append(bits(kmerScore))
				.tab().append(bits(oppositePenalty)).tab();
			text(body, model).tab();
			text(body, anticodon).tab();
			text(body, family).tab().append(partial).tab().append(bits(acceptedRrnaIdentity))
				.tab().append(deferred ? 1 : 0).tab().append(reportStart).tab().append(reportStop)
				.tab().append(proteinEvidence ? 1 : 0).tab().append(ssuEvidence ? 1 : 0).nl();
			if(NcrnaFamilyRegistry.isCobalamin(family)){targets++;}
		}
		@Override
		public void scoreBasis(long id, String family, String model, int nativeModel, int alignedLength,
				float identity, float additive, float multiplier, float assignedScore){
			body.append("B\t").append(id).tab();
			text(body, family).tab();
			text(body, model).tab().append(nativeModel).tab().append(alignedLength).tab().append(bits(identity))
				.tab().append(bits(additive)).tab().append(bits(multiplier)).tab().append(bits(assignedScore)).nl();
			basisRows++;
		}
		@Override
		public void scored(long id, float plus, float minus, int plusLength, int minusLength, long plusPrev, long minusPrev){
			body.append("S\t").append(id).tab().append(bits(plus)).tab().append(bits(minus))
				.tab().append(plusLength).tab().append(minusLength).tab().append(plusPrev).tab().append(minusPrev).nl();
		}
		@Override
		public void end(final long best, final float score, final int selected_){
			require(begun && !ended, "Capture closes the actual native path once");
			selected=selected_;
			body.append("E\t").append(best).tab().append(bits(score)).tab().append(selected).nl();
			ended=true;
		}
		/** The returned buffer is borrowed until the next start; the output writer copies it. */
		public ByteBuilder finish(final ArrayList<Orf> report){
			require(active && ended && report!=null && report.size()==selected && targets==basisRows,
				"Complete native report and one exact basis for every target candidate are required");
			policy.validate(caller);
			output.append(targets>0 ? "GRAPH\t1\t" : "UNCHANGED\t1\t");
			text(output, source).tab().append(read).tab();
			text(output, contig).tab().append(length).tab().append(targets>0 ? count : selected).nl();
			if(targets>0){output.append(body);}
			for(final Orf node:report){output.append("R\t"); node.appendGff(output).nl();}
			output.append("DONE\t").append(selected).nl();
			affected=targets>0;
			active=false;
			return output;
		}
		public GeneCallerPathPolicy policy(){return policy;}
		public boolean affected(){return affected;}
		public int candidates(){return affected ? count : 0;}
		public int targetCandidates(){return targets;}
		@Override
		public void edge(long node, long predecessor, int predecessorStrand, int overlap, float predecessorScore,
				int predecessorLength, float overlapScore, float beforeTransitions, float proposedScore, int proposedLength, boolean replaces){}
		@Override
		public void window(long node, int predecessorList, int lastPosition, int bestPosition, int nominalMinimum, int lastVisited){}
		private final String source;
		private final GeneCaller caller;
		private final ByteBuilder body=new ByteBuilder(), output=new ByteBuilder();
		private GeneCallerPathPolicy policy;
		private String contig;
		private long read;
		private int length, count, selected, targets, basisRows;
		private boolean active, begun, ended, affected;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Saved reader         ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Graph{
		Graph(final String source_, final long read_, final String contig_, final int length_,
				final GeneCallerPathPolicy policy_, final boolean affected, final int count){
			source=source_;
			read=read_;
			contig=contig_;
			length=length_;
			policy=policy_;
			snapshot=affected ? new NcrnaPathSnapshot(source, read, contig) : null;
			headerCount=count;
		}
		public final String source, contig;
		public final long read;
		public final int length;
		public final GeneCallerPathPolicy policy;
		public final NcrnaPathSnapshot snapshot;
		public final ByteBuilder report=new ByteBuilder();
		public final LongList basisNodes=new LongList();
		public final IntList nativeModels=new IntList(), alignedLengths=new IntList();
		public final FloatList identities=new FloatList(), additives=new FloatList(), multipliers=new FloatList();
		public final ArrayList<String> families=new ArrayList<String>(), models=new ArrayList<String>();
		public int reports;
		private final int headerCount;
	}

	/** Owns only parsing state; the supplied ByteFile owner retains stream ownership. */
	public static final class Reader{
		public Reader(final ByteFile input_, final Policies policies_, final int inputs_){
			this(input_, policies_, inputs_, null);
		}
		/** A supplied binding identity requires the bound format; no unbound fallback. */
		public Reader(final ByteFile input_, final Policies policies_, final int inputs_, final String bindingSha80){
			require(input_!=null && policies_!=null && inputs_>0, "Saved graphs require their declared source/policy roster");
			input=input_;
			policies=policies_;
			inputs=inputs_;
			bound=bindingSha80!=null;
			seen=new LongHashSet[inputs];
			final byte[] header=input.nextLine();
			final String expected=bindingSha80==null ? HEADER : NcrnaPathBindings.graphHeader(bindingSha80);
			require(header!=null && arraysEqual(header, expected), "Unsupported or mismatched saved-path binding/header");
		}
		public Graph next(){
			if(complete){return null;}
			final byte[] line=input.nextLine();
			require(line!=null, "Saved paths ended before their complete-file trailer");
			p.set(line);
			if(p.termEquals("COMPLETE", 0)){
				require(p.terms()==5 && p.parseInt(1)==inputs && p.parseLong(2)==contigs
					&& p.parseLong(3)==affected && p.parseLong(4)==reports && input.nextLine()==null,
					"Saved file trailer must reconcile all source/contig/report totals and end the file");
				complete=true;
				return null;
			}
			final boolean graph=p.termEquals("GRAPH", 0);
			require(p.terms()==7 && (graph || p.termEquals("UNCHANGED", 0)) && p.parseInt(1)==1,
				"Expected a versioned complete-contig record");
			final String source=string(p, 2);
			final int index=policies.index(source);
			require(index>=0 && index<inputs, "Saved contig refers to an undeclared source");
			final long read=p.parseLong(3);
			final String contig=string(p, 4);
			final int length=p.parseInt(5), count=p.parseInt(6);
			require(read>=0 && length>=0 && count>=0, "Saved identity, length and candidate counts must be physical");
			if(seen[index]==null){seen[index]=new LongHashSet();}
			require(seen[index].add(read), "A source-qualified read cannot be captured twice");
			final GeneCallerPathPolicy policy=policies.policy(index);
			require(policy!=null, "Saved replay needs the actual input's loaded PGM/policy binding");
			final Graph out=new Graph(source, read, contig, length, policy, graph, count);
			if(graph){readGraph(out);}else{readUnchanged(out, count);}
			contigs++;
			if(graph){affected++;}
			reports+=out.reports;
			return out;
		}
		private void readGraph(final Graph out){
			byte[] row=input.nextLine();
			require(row!=null, "Missing policy before candidate rows");
			p.set(row);
			out.policy.validateRecord(p);
			out.snapshot.begin(out.headerCount==0 ? null : out.contig, out.length, out.headerCount);
			boolean pending=false, ended=false, hasClass=false, hasBasis=false;
			int kind=-1, targets=0, expectedScore=0;
			long current=-1;
			String currentFamily=null, currentModel=null;
			while((row=input.nextLine())!=null){
				p.set(row);
				if(p.termEquals("C", 0)){
					require(!pending && !ended && p.terms()==12, "Candidate rows must precede scoring with complete scalar columns");
					current=p.parseLong(1);
					expectedScore=p.parseInt(8);
					out.snapshot.candidate(current, p.parseInt(2), p.parseInt(3), p.parseInt(4), p.parseInt(5),
						p.parseInt(6), p.parseInt(7), value(expectedScore), value(p.parseInt(9)), flag(p, 10), value(p.parseInt(11)));
					pending=true;
					hasClass=hasBasis=false;
				}else if(p.termEquals("K", 0)){
					require(pending && !hasClass && p.terms()==3 && p.parseLong(1)==current,
						"Each candidate requires one explicit supported subtype");
					kind=p.parseInt(2);
					require(kind==0 || kind==1, "Unknown candidate subtype in saved graph");
					out.snapshot.nodeClass(current, kind==0 ? Orf.class.getName() : NcrnaReportedOrf.class.getName());
					hasClass=true;
				}else if(p.termEquals("M", 0)){
					require(pending && hasClass && p.terms()==21 && p.parseLong(1)==current && flag(p, 16)==(kind==1),
						"Complete metadata must match the immediately preceding candidate and subtype");
					currentModel=nullable(p, 11);
					currentFamily=nullable(p, 13);
					out.snapshot.metadata(current, p.parseInt(2), p.parseInt(3), p.parseInt(4), flag(p, 5), flag(p, 6),
						value(p.parseInt(7)), value(p.parseInt(8)), value(p.parseInt(9)), value(p.parseInt(10)),
						currentModel, nullable(p, 12), currentFamily, p.parseInt(14), value(p.parseInt(15)),
						flag(p, 16), p.parseInt(17), p.parseInt(18), flag(p, 19), flag(p, 20));
					if(NcrnaFamilyRegistry.isCobalamin(currentFamily)){targets++;}
					pending=false;
				}else if(p.termEquals("B", 0)){
					require(!pending && !ended && !hasBasis && p.terms()==10 && p.parseLong(1)==current
						&& NcrnaFamilyRegistry.isCobalamin(currentFamily), "One exact score basis must follow its own target metadata");
					final String family=string(p, 2), model=string(p, 3);
					require(family.equals(currentFamily) && (currentModel==null || model.equals(currentModel))
						&& p.parseInt(4)>=0 && p.parseInt(5)>0 && p.parseInt(9)==expectedScore,
						"Saved model owner/label/span/score must match the native candidate assignment");
					for(int field=6; field<=9; field++){require(Float.isFinite(value(p.parseInt(field))), "Nonfinite saved score basis");}
					if(bound){policies.verifyScoreBasis(family, model, p.parseInt(4));}
					out.basisNodes.add(current);
					out.families.add(family);
					out.models.add(model);
					out.nativeModels.add(p.parseInt(4));
					out.alignedLengths.add(p.parseInt(5));
					out.identities.add(value(p.parseInt(6)));
					out.additives.add(value(p.parseInt(7)));
					out.multipliers.add(value(p.parseInt(8)));
					hasBasis=true;
				}else if(p.termEquals("S", 0)){
					require(!pending && !ended && p.terms()==8, "Scored-state rows must follow complete metadata");
					out.snapshot.scored(p.parseLong(1), value(p.parseInt(2)), value(p.parseInt(3)),
						p.parseInt(4), p.parseInt(5), p.parseLong(6), p.parseLong(7));
				}else if(p.termEquals("E", 0)){
					require(!pending && !ended && p.terms()==4 && out.basisNodes.size==targets && targets>0,
						"Complete target graphs require every exact basis and one native path endpoint");
					out.snapshot.end(p.parseLong(1), value(p.parseInt(2)), p.parseInt(3));
					ended=true;
				}else if(p.termEquals("R", 0)){
					require(ended, "Native report rows follow the complete path oracle");
					report(out, row);
				}else if(p.termEquals("DONE", 0)){
					require(ended && p.terms()==2 && p.parseInt(1)==out.reports, "Closed graph must conserve complete report rows");
					out.snapshot.acceptBaselineRows(out.report, out.reports);
					return;
				}else{throw new IllegalArgumentException("Unknown saved graph record");}
			}
			throw new IllegalArgumentException("Saved candidate graph is truncated");
		}
		private void readUnchanged(final Graph out, final int count){
			for(byte[] row; (row=input.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("R", 0)){report(out, row);}
				else{
					require(p.terms()==2 && p.termEquals("DONE", 0) && p.parseInt(1)==count && out.reports==count,
						"Unchanged contigs retain exactly their original report rows");
					return;
				}
			}
			throw new IllegalArgumentException("Unchanged-contig report is truncated");
		}
		private void report(final Graph out, final byte[] line){
			require(line.length>2 && line[0]=='R' && line[1]=='\t' && p.terms()==10
				&& p.parseInt(4)>0 && p.parseInt(5)>=p.parseInt(4) && p.parseInt(5)<=out.length,
				"Native report framing must preserve a complete physical nine-column GFF row");
			out.report.append(line, 2, line.length-2).nl();
			out.reports++;
		}
		public boolean complete(){return complete;}
		private final ByteFile input;
		private final Policies policies;
		private final int inputs;
		private final boolean bound;
		private final LongHashSet[] seen;
		private final LineParser1 p=new LineParser1('\t');
		private long contigs, affected, reports;
		private boolean complete;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Text primitives        ----------------*/
	/*--------------------------------------------------------------*/

	/** UTF-8 into the reusable byte buffer, without a temporary per-field byte array. */
	public static ByteBuilder text(final ByteBuilder out, final String text){
		if(text==null){return out;}
		for(int i=0; i<text.length(); i++){
			final int c=text.codePointAt(i);
			require(c!=0 && c!='\t' && c!='\n' && c!='\r' && (c<0xd800 || c>0xdfff), "Record identifiers cannot contain delimiters or isolated surrogates");
			if(c<0x80){out.append((byte)c);}
			else if(c<0x800){out.append((byte)(0xc0|(c>>>6))).append((byte)(0x80|(c&63)));}
			else if(c<0x10000){out.append((byte)(0xe0|(c>>>12))).append((byte)(0x80|((c>>>6)&63))).append((byte)(0x80|(c&63)));}
			else{
				out.append((byte)(0xf0|(c>>>18))).append((byte)(0x80|((c>>>12)&63)))
					.append((byte)(0x80|((c>>>6)&63))).append((byte)(0x80|(c&63)));
				i++;
			}
		}
		return out;
	}
	public static String string(final LineParser1 p, final int field){
		final String text=nullable(p, field);
		require(text!=null && !text.isEmpty(), "Required graph identity is empty");
		return text;
	}
	private static String nullable(final LineParser1 p, final int field){
		p.setBounds(field);
		return p.a()==p.b() ? null : new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
	}
	private static boolean flag(final LineParser1 p, final int field){
		final int flag=p.parseInt(field);
		require(flag==0 || flag==1, "Saved boolean flags must be exactly zero or one");
		return flag==1;
	}
	private static int bits(final float value){return Float.floatToRawIntBits(value);}
	private static float value(final int bits){return Float.intBitsToFloat(bits);}
	private static boolean arraysEqual(final byte[] line, final String text){
		if(line.length!=text.length()){return false;}
		for(int i=0; i<line.length; i++){if(line[i]!=text.charAt(i)){return false;}}
		return true;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	public static final String HEADER="#NCRNA_PATH_GRAPHS\t1";
}
