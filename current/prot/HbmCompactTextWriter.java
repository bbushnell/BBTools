package prot;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parser;
import structures.ByteBuilder;

/** Writer-only compact HBM experiment; never changes the input or a runtime model.
 * The accepted dense text supplies all counts. Alphabet order changes storage only.
	 * @author Collei, Keqing
 */
public final class HbmCompactTextWriter {

	public static void main(String[] args){
		String in=null, out=null, order="frequency", filter="none", stats=null, profile=null;
		final String[] expanded=Parser.parseConfig(args);
		if(expanded.length==1 && expanded[0].equalsIgnoreCase("selftest=t")){selfTest(); return;}
		for(String arg:expanded){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Require in=accepted.hbmt.gz out=fresh.txt order=frequency|rows|alphabetical");}
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("in") && in==null){in=value;}
			else if(key.equals("out") && out==null){out=value;}
			else if(key.equals("order")){order=value.toLowerCase(Locale.ROOT);}
			else if(key.equals("filter")){filter=value.toLowerCase(Locale.ROOT);}
			else if(key.equals("stats") && stats==null){stats=value;}
			else if(key.equals("profile") && profile==null){profile=value;}
			else{throw new IllegalArgumentException("Unknown or duplicate option: "+key);}
		}
		if(in==null || out==null || !(order.equals("frequency") || order.equals("rows") || order.equals("alphabetical"))){
			throw new IllegalArgumentException("Require input, fresh output and a supported ordering");
		}
		if(!(filter.equals("none") || filter.equals("singletons") || filter.equals("states") || filter.equals("both"))){
			throw new IllegalArgumentException("filter=none|singletons|states|both is required");
		}
		final HbmCompactTextWriter writer=new HbmCompactTextWriter(order, filter);
		if(profile!=null){writer.annotations=HbmCompactMetadata.fromProfile(profile);}
		writer.process(in, out);
		if(stats!=null){writer.totals.write(stats, filter);}
	}

	private HbmCompactTextWriter(String order_, String filter_){order=order_; filter=filter_;}

	/** Writes owned count graphs through the same family serializer as the dense-input CLI.
	 * Graphs are read synchronously and never modified. Call finish only after all families.
	 * Provenance must describe the actual source; this API does not infer calibration cutoffs.
	 */
	public static final class GraphWriter implements AutoCloseable {
		public GraphWriter(String path, int families, String[] provenance){
			if(families<1 || families>100000 || provenance==null || provenance.length!=16){
				throw new IllegalArgumentException("Graph writer requires1..100000 families and16 provenance pins");
			}
			for(String pin:provenance){DigestSuffix.requireSuffix(pin, "graph provenance");}
			expected=families; output=new HbmCompactTextReader.Output(path);
			preamble(output, "none");
			output.print("#n_families\t"+families+"\n");
			for(int i=0; i<16; i++){output.print("#provenance_"+i+"\t"+provenance[i]+"\n");}
		}
		public void add(String name, AAGraph graph, HbmCompactMetadata metadata){
			if(finished || names.size()>=expected || name==null || name.isEmpty() || !names.add(name)){
				throw new IllegalArgumentException("Duplicate, empty, extra or late graph family");
			}
			final Family family=fromGraph(IdentifierCodec.encode(name), graph);
			if((nodes+=family.rows.size()+(long)graph.ref.length)>HbmBundleFormat.CAP_TOTAL_NODES){
				throw new IllegalArgumentException("Graph bundle exceeds native node cap");
			}
			if(metadata!=null && (metadata.start>=graph.ref.length || metadata.stop>=graph.ref.length)){
				throw new IllegalArgumentException("Metadata endpoint exceeds graph length");
			}
			family.filter(false, false, new Stats()); family.write(output, "frequency", metadata);
		}
		public void finish(){
			if(finished || names.size()!=expected){throw new IllegalStateException("Graph family count differs from declaration");}
			output.finish(); finished=true;
		}
		@Override public void close(){if(output.close()){throw new IllegalStateException("Graph writer output failed");}}
		private final HbmCompactTextReader.Output output;
		private final HashSet<String> names=new HashSet<String>();
		private final int expected;
		private long nodes;
		private boolean finished;
	}

	/** Reject every graph feature the native compact reader cannot reconstruct exactly. */
	private static Family fromGraph(String name, AAGraph graph){
		if(graph==null || graph.ref.length<1 || graph.ref.length>HbmBundleFormat.CAP_L
				|| graph.pivot.length!=graph.ref.length || graph.del.length!=graph.ref.length){
			throw new IllegalArgumentException("Invalid graph dimensions");
		}
		HbmBundleLoader.checkGraphKnobs(graph);
		for(byte residue:graph.pivot){
			if(residue<0 || residue==20 || residue>Blosum62.X_CODE){throw new IllegalArgumentException("Invalid pivot residue");}
		}
		final Family family=new Family(name, AAGraph.decode(graph.pivot).getBytes(java.nio.charset.StandardCharsets.US_ASCII));
		for(int pos=0; pos<graph.ref.length; pos++){
			final AAGraphNode ref=graph.ref[pos], del=graph.del[pos];
			final Row row=graphRow(ref, false, pos, graph.pivot[pos]);
			if(del==null || del.type!=AAGraphNode.DEL || del.rpos!=pos || del.refResidue!=graph.pivot[pos]
					|| del.insEdge!=null || del.count!=null || del.weight!=null || del.countSum<0
					|| del.weightSum!=del.countSum || (long)row.sum+del.countSum>Integer.MAX_VALUE
					|| ((pos==0 || pos==graph.ref.length-1) && del.countSum!=0)){
				throw new IllegalArgumentException("Unrepresentable deletion at "+pos);
			}
			row.deletion=del.countSum; family.rows.add(row);
			int chain=0;
			for(AAGraphNode node=ref.insEdge; node!=null; node=node.insEdge){
				if(pos+1==graph.ref.length || ++chain>HbmBundleFormat.CAP_CHAIN){throw new IllegalArgumentException("Unrepresentable insertion chain at "+pos);}
				family.rows.add(graphRow(node, true, pos+1, Blosum62.X_CODE));
			}
		}
		family.anchor=graph.ref.length-1; return family;
	}

	/** Counts, weights and sums must survive the reader's weight=count reconstruction. */
	private static Row graphRow(AAGraphNode node, boolean insertion, int position, byte residue){
		if(node==null || node.type!=(insertion ? AAGraphNode.INS : AAGraphNode.REF) || node.rpos!=position
				|| node.refResidue!=residue || node.count==null || node.weight==null
				|| node.count.length!=22 || node.weight.length!=22 || node.countSum<1 || node.weightSum!=node.countSum){
			throw new IllegalArgumentException("Unrepresentable graph node at "+position);
		}
		final Row row=new Row(insertion, node.countSum); long sum=0;
		for(int i=0; i<22; i++){
			final int n=node.count[i];
			if(n<0 || node.weight[i]!=n || (i==20 && n!=0)){throw new IllegalArgumentException("Invalid residue count/weight at "+position);}
			row.counts[i]=n; sum+=n;
		}
		if(sum!=row.sum){throw new IllegalArgumentException("Graph histogram sum differs at "+position);}
		return row;
	}

	private static void preamble(HbmCompactTextReader.Output output, String filter){
		output.print("#format\t"+HbmCompactTextReader.FORMAT+"\n#graph_contract\t"+HbmDenseTextBundle.CONTRACT+
			"\n#cutoff_units\t"+HbmCompactTextReader.CUTOFF_UNITS+"\n#coordinates\tzero_based_inclusive\n#sums\toriginal\n#filter\t"+filter+"\n");
		output.print("#columns\t[+ insertion marker]\tsum\tcounts in #alphabet order; trailing zeros omitted\t[-deletion_count]\n");
	}

	/** Buffers only one family so its complete observed alphabet determines every row. */
	private void process(String in, String out){
		final ByteFile input=ByteFile.makeByteFile(in, false);
		final HbmCompactTextReader.Output output=new HbmCompactTextReader.Output(out);
		preamble(output, filter);
		final LineParser1 fields=new LineParser1('\t');
		Family family=null;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				fields.set(line);
				if(fields.terms()<1){throw new IllegalArgumentException("Empty source row");}
				if(fields.termEquals("#format", 0)){header(fields, "hbm_dense_v2"); format++;}
				else if(fields.termEquals("#encoding", 0)){header(fields, "dense"); dense++;}
				else if(fields.termEquals("#coordinate_mode", 0)){header(fields, "explicit"); explicit++;}
				else if(fields.termEquals("#min_count_emitted", 0)){header(fields, "1"); minimum++;}
				else if(fields.termEquals("#n_families", 0)){
					if(fields.terms()!=2 || declared++!=0){throw new IllegalArgumentException("Missing or duplicate family count");}
					expectedFamilies=fields.parseInt(1);
					if(expectedFamilies<1 || expectedFamilies>100000){throw new IllegalArgumentException("Family count outside supported bounds");}
					output.print(new ByteBuilder().append(line).nl());
				}else if(fields.termStartsWith("#provenance_", 0)){
					output.print(new ByteBuilder().append(line).nl());
				}else if(fields.termEquals('f', 0)){
					if(family!=null || fields.terms()!=6 || !names.add(fields.parseString(2))){throw new IllegalArgumentException("Duplicate or unfinished family");}
					if(format!=1 || dense!=1 || explicit!=1 || minimum!=1 || declared!=1){throw new IllegalArgumentException("Source header is not the accepted dense contract");}
					final int length=HbmDenseTextLoader.integer(fields, 3);
					final byte[] consensus=fields.parseByteArray(5);
					if(length<1 || consensus.length!=length){throw new IllegalArgumentException("Consensus length mismatch");}
					family=new Family(fields.parseString(2), consensus);
				}else if(fields.termEquals('e', 0)){
					if(family==null || fields.terms()!=2){throw new IllegalArgumentException("Expected v2 family terminator");}
					DigestSuffix.requireSuffix(fields.parseString(1), "source family checksum");
					family.filter(filter.equals("singletons") || filter.equals("both"),
						filter.equals("states") || filter.equals("both"), totals);
					final HbmCompactMetadata info=annotations==null ? null : annotations.remove(IdentifierCodec.decode(family.name));
					if(annotations!=null && info==null){throw new IllegalArgumentException("Profile missing family "+family.name);}
					family.write(output, order, info);
					refRows+=family.consensus.length; insRows+=family.retainedRows-family.consensus.length;
					family=null;
				}else if(line[0]!='#' && !fields.termEquals('z', 0)){
					if(family==null){throw new IllegalArgumentException("Data row outside family");}
					family.add(fields);
				}
			}
			if(family!=null || names.size()!=expectedFamilies || (annotations!=null && !annotations.isEmpty())){
				throw new IllegalArgumentException("Incomplete or mismatched source/profile families: "+names.size());
			}
			output.finish();
		}finally{
			final boolean readError=input.close(), writeError=output.close();
			if(readError || writeError){throw new IllegalStateException("Compact writer I/O failed");}
		}
		System.err.println("HBM_COMPACT_WRITE_PASS families="+names.size()+" reference_rows="+refRows+" insertion_rows="+insRows+" order="+order+" filter="+filter);
	}

	private static void header(LineParser1 fields, String value){
		if(fields.terms()!=2 || !fields.termEquals(value, 1)){throw new IllegalArgumentException("Unexpected header; expected "+value);}
	}

	/** Native residue slots remain unchanged; zero-frequency symbols need no column. */
	private static final class Family {
		Family(String name_, byte[] consensus_){name=name_; consensus=consensus_;}

		void add(LineParser1 fields){
			if(fields.terms()<3){throw new IllegalArgumentException("Short HBM data row");}
			final int position=HbmDenseTextLoader.integer(fields, 1);
			if(fields.termEquals('r', 0)){
				if(position!=anchor+1 || position>=consensus.length){throw new IllegalArgumentException("Reference order changed");}
				anchor=position; chain=0; deletionSeen=false;
				reference=counts(fields, 2, false); rows.add(reference);
			}else if(fields.termEquals('i', 0)){
				if(reference==null || position!=anchor || deletionSeen || (long)position+1>=consensus.length
						|| HbmDenseTextLoader.integer(fields, 2)!=chain++){
					throw new IllegalArgumentException("Insertion order or boundary is not representable");
				}
				rows.add(counts(fields, 3, true));
			}else if(fields.termEquals('d', 0)){
				if(reference==null || position!=anchor || deletionSeen || fields.terms()!=3 || position==0 || position==consensus.length-1){
					throw new IllegalArgumentException("Deletion order or boundary is not representable");
				}
				reference.deletion=HbmDenseTextLoader.integer(fields, 2); deletionSeen=true;
				if(reference.deletion<1){throw new IllegalArgumentException("Empty DEL row must be omitted");}
			}else if(fields.termEquals("di", 0)){
				throw new IllegalArgumentException("DEL-attached insertion requires a format decision; none may be discarded");
			}else{throw new IllegalArgumentException("Unknown HBM row type");}
		}

		/** Reuses the native canonical A48 decoder and retains the original sum exactly. */
		private Row counts(LineParser1 fields, int mode, boolean insertion){
			if(fields.terms()!=mode+24 || !fields.termEquals('d', mode)){throw new IllegalArgumentException("Expected22 dense residue counts");}
			final Row row=new Row(insertion, HbmDenseTextLoader.integer(fields, mode+1));
			long observed=0;
			for(int i=0; i<22; i++){
				final int count=row.counts[i]=HbmDenseTextLoader.integer(fields, mode+2+i);
				observed+=count;
			}
			if(row.sum<1 || observed!=row.sum){throw new IllegalArgumentException("Original countSum and histogram differ");}
			return row;
		}

		/** Uses original positional depths, then rebuilds alphabet frequencies after filtering.
		 * INS follows the earlier HbmRareStateProbe.pruneChain rule: denominator at
		 * anchor+1 and removal of the complete rare suffix, never shifting nodes.
		 * Singleton filtering retains the original countSum even when counts fall.
		 */
		void filter(boolean singletons, boolean states, Stats stats){
			final long[] depth=new long[consensus.length]; int position=-1;
			for(Row row:rows){
				stats.inputRows++;
				for(int count:row.counts){if(count>0){stats.inputResidueCells++;}}
				if(!row.insertion){depth[++position]=(long)row.sum+row.deletion;}
				if(row.deletion>0){stats.inputDeletionTokens++;}
			}
			if(position!=depth.length-1){throw new IllegalArgumentException("Reference row count differs from consensus");}
			position=-1; boolean suffix=false;
			for(Row row:rows){
				if(!row.insertion){position++; suffix=false;}
				if(states && row.insertion){
					if(position+1>=depth.length){throw new IllegalArgumentException("Insertion has no following reference depth");}
					suffix|=rare(row.sum, depth[position+1]);
					if(suffix){
						row.removed=true; stats.removedRows++;
						for(int count:row.counts){if(count>0){stats.droppedRowCells++;}}
						continue;
					}
				}
				if(states && row.deletion>0 && rare(row.deletion, depth[position])){
					row.deletion=0; stats.removedDeletionTokens++;
				}
				retainedRows++; stats.outputRows++;
				for(int i=0; i<row.counts.length; i++){
					if(singletons && row.sum>=100 && row.counts[i]==1){row.counts[i]=0; stats.singletonCells++;}
					mass[i]+=row.counts[i];
					if(row.counts[i]>0){nonzeroRows[i]++; stats.outputResidueCells++;}
				}
			}
		}

		void write(HbmCompactTextReader.Output output, String order, HbmCompactMetadata metadata){
			if(anchor!=consensus.length-1){throw new IllegalArgumentException("Missing reference positions");}
			final int[] alphabet=new int[22]; int size=0;
			for(int i=0; i<22; i++){if(mass[i]>0){alphabet[size++]=i;}}
			for(int i=1; i<size; i++){
				final int symbol=alphabet[i]; int j=i;
				while(j>0 && before(symbol, alphabet[j-1], order)){alphabet[j]=alphabet[j-1]; j--;}
				alphabet[j]=symbol;
			}
			assert(size>0) : "Native pivot scaffold gives each HBM at least one observed residue (AAGraph constructor)";
			final ByteBuilder buffer=new ByteBuilder(4096);
			buffer.append("#name\t").append(name).nl().append("#alphabet\t");
			for(int i=0; i<size; i++){buffer.append(symbol(alphabet[i]));}
			buffer.nl().append("#consensus\t").append(consensus).nl();
			if(metadata!=null){metadata.append(buffer);}
			buffer.append("#rows\t").append(retainedRows).nl();
			output.print(buffer); buffer.clear();
			for(Row row:rows){
				if(row.removed){continue;}
				if(row.insertion){buffer.append('+').tab();}
				buffer.appendA48(row.sum);
				int last=size-1;
				while(last>=0 && row.counts[alphabet[last]]==0){last--;}
				for(int i=0; i<=last; i++){buffer.tab().appendA48(row.counts[alphabet[i]]);}
				if(row.deletion>0){buffer.tab().append('-').appendA48(row.deletion);}
				buffer.nl(); output.print(buffer); buffer.clear();
			}
		}

		private boolean before(int a, int b, String order){
			final long x=order.equals("frequency") ? mass[a] : order.equals("rows") ? nonzeroRows[a] : 0;
			final long y=order.equals("frequency") ? mass[b] : order.equals("rows") ? nonzeroRows[b] : 0;
			return x>y || (x==y && symbol(a)<symbol(b));
		}

		final String name;
		final byte[] consensus;
		final ArrayList<Row> rows=new ArrayList<Row>();
		final long[] mass=new long[22], nonzeroRows=new long[22];
		int anchor=-1, chain, retainedRows;
		Row reference;
		boolean deletionSeen;
	}

	private static char symbol(int code){return code==Blosum62.X_CODE ? 'X' : (char)AminoAcid.numberToAcid[code];}

	/** Exact strict1% cutoff, with long arithmetic; equality is retained. */
	private static boolean rare(int count, long depth){return count>0 && 100L*count<depth;}

	/** Boundary fixture checks the original denominator, suffix topology and unchanged sums. */
	private static void selfTest(){
		final Family family=new Family("fixture", new byte[]{'A', 'A', 'A'});
		final Row first=fixtureRow(false, 100), ins1=fixtureRow(true, 2), ins2=fixtureRow(true, 1);
		final Row middle=fixtureRow(false, 199), exact=fixtureRow(true, 1), last=fixtureRow(false, 99);
		middle.deletion=2;
		family.rows.add(first); family.rows.add(ins1); family.rows.add(ins2);
		family.rows.add(middle); family.rows.add(exact); family.rows.add(last);
		final Stats stats=new Stats(); family.filter(true, true, stats);
		if(!ins1.removed || !ins2.removed || exact.removed || middle.deletion!=0
				|| first.sum!=100 || first.counts[1]!=0 || middle.sum!=199 || middle.counts[1]!=0
				|| last.counts[1]!=1 || stats.singletonCells!=2 || stats.removedRows!=2
				|| stats.droppedRowCells!=2 || stats.removedDeletionTokens!=1
				|| family.retainedRows!=4 || family.mass[1]!=1 || family.nonzeroRows[1]!=1
				|| rare(1, 100) || !rare(1, 101) || rare(Integer.MAX_VALUE, Integer.MAX_VALUE)){
			throw new AssertionError("Compact filters changed cutoff, original sums, next-position depth or suffix topology");
		}
		stats.validate();
		System.err.println("HBM_COMPACT_FILTER_FIXTURE_PASS");
	}

	private static Row fixtureRow(boolean insertion, int sum){
		final Row row=new Row(insertion, sum);
		row.counts[0]=sum;
		if(!insertion){row.counts[0]--; row.counts[1]=1;}
		return row;
	}

	/** Disjoint removal categories make the combined arm's accounting additive. */
	private static final class Stats {
		void validate(){
			if(inputRows!=outputRows+removedRows || inputResidueCells!=outputResidueCells+singletonCells+droppedRowCells){
				throw new IllegalStateException("Writer removal accounting failed conservation");
			}
		}
		void write(String path, String filter){
			validate();
			final ByteBuilder text=new ByteBuilder("filter\tinput_rows\toutput_rows\tremoved_insertion_rows\tinput_residue_cells\toutput_residue_cells\tsingleton_cells_removed\tdropped_row_residue_cells\tinput_deletion_tokens\tremoved_deletion_tokens\tall_nonzero_cells_removed\trule\n");
			text.append(filter).tab().append(inputRows).tab().append(outputRows).tab().append(removedRows).tab()
				.append(inputResidueCells).tab().append(outputResidueCells).tab().append(singletonCells).tab()
				.append(droppedRowCells).tab().append(inputDeletionTokens).tab().append(removedDeletionTokens).tab()
				.append(singletonCells+droppedRowCells+removedDeletionTokens).tab()
				.append("A:count==1 and original sum>=100; B:100*stateCount<original REF+DEL depth, INS at anchor+1 and suffix removal, DEL at anchor; C:B then A; alphabet after filtering; sums unchanged").nl();
			final ByteStreamWriter out=new ByteStreamWriter(path, false, false, false); out.start();
			try{out.print(text);}finally{if(out.poisonAndWait()){throw new IllegalStateException("Filter statistics write failed");}}
		}
		long inputRows, outputRows, removedRows, inputResidueCells, outputResidueCells;
		long singletonCells, droppedRowCells, inputDeletionTokens, removedDeletionTokens;
	}

	private static final class Row {
		Row(boolean insertion_, int sum_){insertion=insertion_; sum=sum_;}
		final boolean insertion;
		final int sum;
		final int[] counts=new int[22];
		int deletion;
		boolean removed;
	}

	private final String order, filter;
	private final Stats totals=new Stats();
	private HashMap<String,HbmCompactMetadata> annotations;
	private final HashSet<String> names=new HashSet<String>();
	private int format, dense, explicit, minimum, declared, expectedFamilies;
	private long refRows, insRows;
}
