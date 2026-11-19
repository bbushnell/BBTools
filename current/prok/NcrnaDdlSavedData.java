package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Consumer;

import fileIO.FileFormat;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.SequenceInput;
import stream.Read;
import stream.StreamerFactory;
import shared.Tools;
import structures.IntHashMap;
import structures.ListNum;

/** Reads the accepted selected inventory and saved ranking cache, never genomes.
 * Original unique attribution is borrowed from links, not recomputed after
 * target/input selection. One owner binds metadata and retains window bytes
 * before workers begin; workers subsequently borrow all state read-only.
 * @author Brian, Raiden
 */
public final class NcrnaDdlSavedData{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static NcrnaDdlSavedData load(final String manifest, final String inventory,
			final int models, final int expectedGroups){
		require(manifest!=null && inventory!=null && models>0 && expectedGroups>0,
			"Saved windows require accepted panels/inventory and complete model/group counts");
		final NcrnaDdlSavedData data=new NcrnaDdlSavedData(models);
		data.panels(manifest);
		data.inventory(inventory, expectedGroups);
		for(final Panel panel:data.panels){
			data.captureInputs(panel);
			data.parents(panel);
			data.loci(panel);
			data.windows(panel);
			data.ranks(panel);
			data.links(panel);
			data.sequences(panel);
		}
		return data;
	}

	private NcrnaDdlSavedData(final int models_){models=models_;}

	/** Original saved order; no new window selection or scoring is performed. */
	public int windowCount(){return windows.size();}
	/** Borrowed read-only sequence bytes for one original saved window. */
	public byte[] windowBases(final int index){return windows.get(index).bases;}

	/** Optional hybrid prepass: reads the original model-complete native counts,
	 * never runs an index or verifier. Full model headers bind each original ordinal
	 * to the unchanged reference FASTA, including its retained description. */
	public void readNativeCounts(final String[] modelNames){
		readNativeCounts(modelNames, null, false);
	}

	/** Optional GC observations use each accepted model/window row once, before
	 * parent links can duplicate it. Statistics retain count,mean,M2,min,max,missing.
	 * Accepted bits are requested only by the GC judging study. */
	public void readNativeCounts(final String[] modelNames, final double[][][] gc, final boolean keepAccepted){
		require(!nativeStarted && modelNames!=null && modelNames.length==models,
			"Read native counts once against the complete original model roster before starting workers");
		nativeStarted=true;
		require(gc==null || gc.length==3 && gc[0].length==models && gc[1].length==models && gc[2].length==models,
			"GC descriptions require the complete member,padded-window,and accepted-span model arrays");
		if(gc!=null || keepAccepted){
			for(final Window window:windows){window.gc=Tools.calcGC(window.bases);}
		}
		final int[] bounds=new int[2];
		final byte[][] names=new byte[models][];
		for(int model=0; model<models; model++){
			require(modelNames[model]!=null && !modelNames[model].isEmpty(), "Cached native ordinals require their full original model headers");
			names[model]=modelNames[model].getBytes(StandardCharsets.UTF_8);
		}
		for(final Panel panel:panels){
			try(Input in=new Input(panel.prefix+".oracle.tsv.gz")){
				header(in, ORACLE_HEADER);
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					parse(row, 15);
					final Window window=event(panel, p);
					if(window==null){continue;}
					final int model=p.parseInt(5), count=p.parseInt(7);
					if(model<0 || model>=models || model!=window.nativeModels || !same(p, 6, names[model])){
						throw new IllegalArgumentException("Every cached window needs each native model once, in original ordinal/full-header order: panel="
							+panel.name+", event="+window.key+", expected_model="+window.nativeModels+", actual_model="+model+", header="+text(p, 6));
					}
					if(count<0 || count>65535){
						throw new IllegalArgumentException("Cached native counts must satisfy TrnaKmerIndex.MAX_COUNT=65535, preserving exact float ranking: panel="
							+panel.name+", event="+window.key+", model="+model+", count="+count);
					}
					if(window.nativeCounts==null){window.nativeCounts=new char[models];}
					window.nativeCounts[model]=(char)count;
					window.nativeModels++;
					if(bool(p, 9)){
						final int start=p.parseInt(10), stop=p.parseInt(11);
						require(start>0 && stop>=start, "Accepted cached model spans are one-based inclusive original-forward coordinates");
						window.acceptedModels++;
						window.acceptedStart=Math.min(window.acceptedStart, start);
						window.acceptedStop=Math.max(window.acceptedStop, stop);
						if(keepAccepted){
							if(window.acceptedModelBits==null){window.acceptedModelBits=new long[(models+63)/64];}
							window.acceptedModelBits[model>>>6]|=1L<<(model&63);
						}
						if(gc!=null){
							addGc(gc[1][model], window.gc);
							BuildNcrnaDdlLayoutCounts.paddingBounds(window.start, window.stop, window.strand, start, stop, bounds);
							addGc(gc[2][model], Tools.calcGC(window.bases, bounds[0], bounds[1]-1));
						}
					}
				}
			}
		}
		for(final Window window:windows){
			if(window.nativeCounts==null || window.nativeModels!=models){
				throw new IllegalArgumentException("Selected windows cannot lose zero-support models while reconstructing the original native order: panel="
					+window.input.panel.name+", event="+window.key+", found="+window.nativeModels+", expected="+models);
			}
		}
	}

	/** Primitive Welford accumulators; no sequence or observation objects. */
	public static double[][] gcStatistics(final int models){
		require(models>0, "GC statistics require an explicit nonempty model roster");
		final double[][] result=new double[models][6];
		for(final double[] row:result){row[3]=Double.POSITIVE_INFINITY; row[4]=Double.NEGATIVE_INFINITY;}
		return result;
	}
	public static void addGc(final double[] stats, final double value){
		assert(stats!=null && stats.length==6) : "Welford GC rows preserve count,mean,M2,min,max,and missing count";
		if(Double.isNaN(value)){stats[5]++; return;}
		require(Double.isFinite(value) && value>=0 && value<=1, "GC fractions use defined A/C/G/T bases only");
		final double delta=value-stats[1];
		stats[0]++;
		stats[1]+=delta/stats[0];
		stats[2]+=delta*(value-stats[1]);
		stats[3]=Math.min(stats[3], value);
		stats[4]=Math.max(stats[4], value);
	}
	public static double gcSd(final double[] stats){
		assert(stats!=null && stats.length==6 && stats[0]>=0) : "GC SD uses the same count and centered sum as its Welford row";
		return stats[0]==0 ? Double.NaN : stats[0]==1 ? 0 : Math.sqrt(Math.max(0, stats[2]/(stats[0]-1)));
	}

	/*--------------------------------------------------------------*/
	/*----------------         Input Selection      ----------------*/
	/*--------------------------------------------------------------*/

	private void panels(final String path){
		try(Input in=new Input(path)){
			header(in, NcrnaSavedRankSummary.PANEL_HEADER);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 11);
				final String name=text(p, 0), prefix=text(p, 1), labels=text(p, 2);
				require(!name.isEmpty() && prefix.startsWith("/") && labels.startsWith("/") && byName.get(name)==null,
					"Each saved panel requires a unique name and absolute capture/label paths");
				final Panel panel=new Panel(name, prefix, labels, p.parseInt(5), p.parseInt(6),
					p.parseInt(7), p.parseInt(8), p.parseInt(9), p.parseInt(10));
				panels.add(panel);
				byName.put(name, panel);
			}
		}
		require(!panels.isEmpty(), "A saved experiment needs at least one declared panel");
	}

	private void inventory(final String path, final int expectedGroups){
		final IntHashMap groups=new IntHashMap();
		final ObjectSet<String> ids=new ObjectSet<String>(String.class), paths=new ObjectSet<String>(String.class);
		try(Input in=new Input(path)){
			header(in, INVENTORY_HEADER);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 11);
				final Panel panel=byName.get(text(p, 0));
				require(panel!=null, "Selected inventory cannot introduce an undeclared capture panel");
				final int slot=p.parseInt(1), group=p.parseInt(4);
				final String id=text(p, 2), fasta=text(p, 3);
				require(slot>=0 && slot<panel.inputs.length && panel.inputs[slot]==null && group>=0
					&& !id.isEmpty() && fasta.startsWith("/") && !groups.contains(group) && ids.add(id) && paths.add(fasta),
					"Selected groups, original IDs/paths and panel slots must be unique and resolved");
				groups.put(group, 1);
				final Selected input=new Selected(selected.size(), panel, slot, id, fasta, group, p);
				panel.inputs[slot]=input;
				panel.selected++;
				selected.add(input);
			}
		}
		require(groups.size()==expectedGroups, "Keep every selected group, including zero-target inputs; actual="+groups.size());
		for(final Panel panel:panels){
			require(panel.selected==panel.expectedSelected, "Selected inventory must conserve each declared panel subset: "+panel.name);
		}
	}

	private void captureInputs(final Panel panel){
		try(Input in=new Input(panel.prefix+".inputs.tsv")){
			header(in, NcrnaDdlTrainingOverlap.SAMPLE_HEADER);
			int count=0;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 13);
				require(count<panel.inputs.length && p.parseInt(0)==count, "Original capture inputs require complete dense slots");
				final Selected input=panel.inputs[count++];
				if(input==null){continue;}
				require(same(p, 1, input.idBytes) && same(p, 2, input.pathBytes) && p.parseLong(5)>=0,
					"Selected input must bind its exact original capture ID/path and window count");
				input.expectedWindows=p.parseLong(5);
			}
			require(count==panel.inputs.length, "The original capture census cannot omit selected zero-target inputs");
		}
	}

	private void parents(final Panel panel){
		final LineParser1 label=new LineParser1('\t');
		try(Input in=new Input(panel.prefix+".parents.tsv.gz"); Input labels=new Input(panel.labels)){
			header(in, NcrnaSavedRankSummary.PARENT_HEADER);
			header(labels, NcrnaSavedRankSummary.STRATA_HEADER);
			int targets=0;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 14);
				if(!bool(p, 9)){continue;}
				targets++;
				final byte[] labelRow=labels.reader.nextLine();
				require(labelRow!=null, "Saved labels must retain every original target before subset selection");
				label.set(labelRow);
				require(label.terms()==4 && label.parseInt(0)==p.parseInt(0) && sameFields(label, 1, p, 1)
					&& sameFields(label, 2, p, 3), "Saved seen-status labels must retain original parent order and identity");
				final int slot=p.parseInt(0);
				require(slot>=0 && slot<panel.inputs.length, "A target parent must belong to its original capture panel");
				final Selected input=panel.inputs[slot];
				if(input==null){continue;}
				require(!bool(p, 8) && same(p, 1, input.idBytes) && same(p, 2, input.pathBytes),
					"Selected targets cannot promote masks or exchange input identities");
				final String id=text(p, 3);
				require(panel.parents.get(id)==null, "Each selected target needs one original parent record");
				final Parent parent=new Parent(parents.size(), input, id, text(p, 4), p.parseInt(5), p.parseInt(6), text(p, 7), stratum(label));
				panel.parents.put(id, parent);
				parents.add(parent);
				input.parents++;
			}
			require(targets==panel.expectedTargets && labels.reader.nextLine()==null, "All original target labels must be conserved");
		}
		for(final Selected input:panel.inputs){
			if(input!=null){require(input.parents==input.expectedTargets, "Selected target inventory must include every target and zero-target input");}
		}
	}

	private void loci(final Panel panel){
		final LineParser1 first=new LineParser1('\t');
		int count=0, common=0, zero=0, gaps=0;
		try(Input in=new Input(panel.prefix+".loci.tsv.gz")){
			header(in, NcrnaSavedRankSummary.LOCUS_HEADER);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				first.set(row);
				require(first.terms()==12 && first.parseInt(3)==0, "Every original locus block begins with native method0");
				final Parent parent=panel.parents.get(text(first, 2));
				if(parent!=null){
					require(!parent.locusRead && first.parseInt(0)==parent.input.slot && same(first, 1, parent.input.idBytes), "Saved minima cannot duplicate or exchange selected parent bindings");
					parent.locusRead=true;
					parent.nativeRank=first.parseInt(8);
					parent.windows=first.parseInt(4);
					parent.compatible=first.parseInt(5);
					parent.ambiguous=first.parseInt(6);
					parent.unique=first.parseInt(7);
					parent.reason=text(first, 9);
					if(parent.nativeRank>0){parent.input.common++;}
					else if(parent.reason.equals("ZERO_NATIVE_INDEX_SUPPORT")){parent.input.zero++;}
					else if(parent.reason.equals("NO_COMPATIBLE_MODEL")){parent.input.gaps++;}
					else if(parent.reason.equals("NO_CAPTURED_WINDOW")){parent.input.noWindow++;}
					else if(parent.reason.equals("AMBIGUOUS_SPANS_ONLY")){parent.input.ambiguous++;}
				}
				if(first.parseInt(8)>0){common++;}
				if(first.termEquals("ZERO_NATIVE_INDEX_SUPPORT", 9)){zero++;}
				if(first.termEquals("NO_COMPATIBLE_MODEL", 9)){gaps++;}
				for(int method=0; method<METHODS; method++){
					if(method>0){row=in.reader.nextLine();}
					require(row!=null, "Saved locus blocks cannot end before all original methods");
					parse(row, 12);
					require(p.parseInt(0)==first.parseInt(0) && sameFields(p, 1, first, 1) && sameFields(p, 2, first, 2)
						&& p.parseInt(3)==method, "Saved locus minima must keep each complete original method block");
					for(int field=4; field<=7; field++){
						require(p.parseInt(field)==first.parseInt(field), "A saved method cannot change a locus's oracle disposition");
					}
					final int windows=p.parseInt(4), compatible=p.parseInt(5), unique=p.parseInt(7), rank=p.parseInt(8);
					require(windows>=0 && compatible>=0 && compatible<=windows && p.parseInt(6)>=0 && unique>=0
						&& (long)p.parseInt(6)+unique==compatible, "Compatible locus windows partition into unique and ambiguous events");
					final String reason=windows==0 ? "NO_CAPTURED_WINDOW" : disposition(compatible, unique, rank);
					require(p.termEquals(reason, 9) && (reason.equals("RANKED") ? rank>0 && rank<=models : rank==-1),
						"Saved locus rank and reason must agree with RankNcrnaWindows.reason");
					require(method==0 || unique==0 || rank>0, "Every sketch orders all models even when its scores are zero");
					if(parent!=null && CONTROL_INDEX[method]>=0){parent.controlMin[CONTROL_INDEX[method]]=p.parseInt(8);}
				}
				count++;
			}
		}
		require(count==panel.expectedTargets && common==panel.expectedCommon && zero==panel.expectedZero && gaps==panel.expectedGaps,
			"Source minima must conserve the accepted complete target/common/zero/gap counts before selection");
		for(final Parent parent:panel.parents.values()){
			if(parent!=null){require(parent.locusRead, "Every selected parent requires its original saved minima");}
		}
		for(final Selected input:panel.inputs){
			if(input!=null){
				require(input.common==input.expectedCommon && input.zero==input.expectedZero && input.gaps==input.expectedGaps
					&& input.noWindow==input.expectedNoWindow && input.ambiguous==input.expectedAmbiguous,
					"Selected common/zero/gap/no-window/ambiguous counts must match the accepted inventory");
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Window Evidence       ----------------*/
	/*--------------------------------------------------------------*/

	private void windows(final Panel panel){
		try(Input in=new Input(panel.prefix+".windows.tsv.gz")){
			header(in, NcrnaDdlLocusOverlap.WINDOW_HEADER);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 15);
				final Selected input=selectedInput(panel, p);
				if(input==null){continue;}
				final Window window=new Window(input, p.parseLong(2), p.parseLong(3), text(p, 4), p.parseInt(6),
					p.parseInt(10), p.parseInt(11), p.parseInt(12), p.parseInt(14));
				require(panel.windows.get(window.key)==null, "Saved window IDs must be unique within their panel");
				panel.windows.put(window.key, window);
				windows.add(window);
				input.windows++;
			}
		}
		for(final Selected input:panel.inputs){
			if(input!=null){require(input.windows==input.expectedWindows, "Selected window metadata must match the original capture accounting");}
		}
	}

	private void ranks(final Panel panel){
		try(Input in=new Input(panel.prefix+".ranks.tsv.gz")){
			header(in, RANK_HEADER);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 16);
				final Window window=event(panel, p);
				if(window==null){continue;}
				final Edge edge=edge(panel, window, p, 5, true);
				final int method=p.parseInt(6), compatible=p.parseInt(9), unique=p.parseInt(10), rank=p.parseInt(11);
				require(method>=0 && method<METHODS && p.parseInt(7)==window.length
					&& p.parseInt(8)==edge.parent.stop-edge.parent.start+1 && compatible>=0 && compatible<=models
					&& unique>=0 && unique<=compatible, "Saved event ranks must retain their exact window/parent dimensions and oracle counts");
				final int word=method>>>6;
				final long bit=1L<<(method&63);
				require((edge.methods[word]&bit)==0, "A saved event/parent may contain each original method only once");
				edge.methods[word]|=bit;
				edge.methodCount++;
				if(method==0){
					edge.compatible=compatible;
					edge.unique=unique;
					edge.nativeRank=rank;
				}else{require(compatible==edge.compatible && unique==edge.unique, "A ranking method cannot change its cached oracle disposition");}
				final int tieFirst=p.parseInt(12), tieLast=p.parseInt(13), model=p.parseInt(14);
				require(rank==-1 ? tieFirst==-1 && tieLast==-1 && model==-1
					: rank>0 && rank<=models && tieFirst>0 && tieFirst<=rank && tieLast>=rank && tieLast<=models && model>=0 && model<models,
					"Every saved rank must retain its complete tie span and selected original model");
				require(unique>0 ? method==0 || rank>0 : rank==-1, "Only uniquely attributed oracle spans may enter an event ranking");
				require(p.termEquals(disposition(compatible, unique, rank), 15), "Saved event rank/reason must agree with RankNcrnaWindows.rankRow");
				if(method==0){
					edge.nativeTuple[0]=rank;
					edge.nativeTuple[1]=tieFirst;
					edge.nativeTuple[2]=tieLast;
					edge.nativeTuple[3]=model;
				}
				final int control=CONTROL_INDEX[method];
				if(control>=0){
					final int offset=control*4;
					for(int i=0; i<4; i++){edge.controls[offset+i]=p.parseInt(11+i);}
				}
			}
		}
		for(final Window window:windows){
			if(window.input.panel!=panel){continue;}
			for(final Edge edge:window.edges){
				require(edge.methodCount==METHODS, "Every included event/parent requires all original method rows, including empty-oracle cases");
				final Parent parent=edge.parent;
				parent.events++;
				if(edge.compatible>0){parent.compatibleEvents++;}
				if(edge.unique>0){parent.uniqueEvents++;}
				else if(edge.compatible>0){parent.ambiguousEvents++;}
				parent.observedNative=minimum(parent.observedNative, edge.nativeRank);
				for(int c=0; c<CONTROLS; c++){parent.observedControls[c]=minimum(parent.observedControls[c], edge.controls[c*4]);}
			}
		}
		for(final Parent parent:panel.parents.values()){
			if(parent==null){continue;}
			require(parent.events==parent.windows && parent.compatibleEvents==parent.compatible
				&& parent.uniqueEvents==parent.unique && parent.ambiguousEvents==parent.ambiguous
				&& parent.observedNative==parent.nativeRank && Arrays.equals(parent.observedControls, parent.controlMin),
				"Saved event evidence must recombine to every selected original locus minimum and disposition");
		}
	}

	private void links(final Panel panel){
		try(Input in=new Input(panel.prefix+".links.tsv.gz")){
			header(in, LINK_HEADER);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 10);
				final Window window=event(panel, p);
				if(window==null || !bool(p, 7)){continue;}
				final Edge edge=edge(panel, window, p, 6, false);
				final int model=p.parseInt(5), links=p.parseInt(8);
				final boolean unique=bool(p, 9);
				require(model>=0 && model<models && links>0 && unique==(links==1),
					"A target's cached unique link must retain its original all-scored-parent attribution");
				final int word=model>>>6;
				final long bit=1L<<(model&63);
				require((edge.seenModels[word]&bit)==0, "A compatible model cannot vote twice for one saved event/parent");
				edge.seenModels[word]|=bit;
				edge.links++;
				if(unique){
					edge.uniqueModels[word]|=bit;
					edge.uniqueLinks++;
				}
			}
		}
		for(final Window window:windows){
			if(window.input.panel!=panel){continue;}
			for(final Edge edge:window.edges){
				require(edge.links==edge.compatible && edge.uniqueLinks==edge.unique,
					"Cached links must conserve every compatible/unique model declared by the saved rank rows");
			}
		}
	}

	private void sequences(final Panel panel){
		sequences(panel.prefix+".windows.fa.gz", read->{
			final int split=read.id.indexOf('_');
			require(split>0, "Captured FASTA names require an original slot_event identifier");
			final int slot=Integer.parseInt(read.id.substring(0, split));
			require(slot>=0 && slot<panel.inputs.length, "Captured FASTA slot must belong to its original panel");
			if(panel.inputs[slot]==null){return;}
			final Window window=panel.windows.get(read.id);
			require(window!=null && window.bases==null && read.length()==window.length,
				"Selected saved windows must have one exact-length FASTA record each");
			window.bases=read.bases;
		});
		for(final Window window:windows){
			if(window.input.panel==panel){require(window.bases!=null, "Selected window metadata cannot silently lose its sequence");}
		}
	}

	private Selected selectedInput(final Panel panel, final LineParser1 p){
		final int slot=p.parseInt(0);
		require(slot>=0 && slot<panel.inputs.length, "Saved metadata must retain its original input-slot namespace");
		final Selected input=panel.inputs[slot];
		if(input!=null){require(same(p, 1, input.idBytes), "Saved metadata cannot exchange selected original input IDs");}
		return input;
	}
	private Window event(final Panel panel, final LineParser1 p){
		final Selected input=selectedInput(panel, p);
		if(input==null){return null;}
		final long event=p.parseLong(2);
		if(cachedWindow==null || cachedWindow.input!=input || cachedWindow.event!=event){
			cachedWindow=panel.windows.get(input.slot+"_"+event);
		}
		require(cachedWindow!=null && cachedWindow.ordinal==p.parseLong(3) && same(p, 4, cachedWindow.seqidBytes),
			"Saved event metadata must agree on panel, original slot, event, contig ordinal and seqid");
		return cachedWindow;
	}
	private Edge edge(final Panel panel, final Window window, final LineParser1 p, final int field, final boolean create){
		if(cachedEdge!=null && cachedEdge.window==window && same(p, field, cachedEdge.parent.idBytes)){return cachedEdge;}
		final Parent parent=panel.parents.get(text(p, field));
		require(parent!=null && parent.input==window.input && Arrays.equals(parent.seqidBytes, window.seqidBytes)
			&& parent.strand==window.strand && parent.start-1<=window.stop && parent.stop-1>=window.start,
			"Saved target links/ranks must bind the selected parent to the same input, contig, strand and overlapping window");
		for(final Edge edge:window.edges){
			if(edge.parent==parent){
				cachedEdge=edge;
				return edge;
			}
		}
		require(create && p.parseInt(6)==0, "A compatible link requires its original native-first event-rank block");
		cachedEdge=new Edge(window, parent, models);
		window.edges.add(cachedEdge);
		return cachedEdge;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Utilities          ----------------*/
	/*--------------------------------------------------------------*/

	/** Streams bound reference/window FASTA through native I/O. The consumer
	 * uses borrowed records synchronously before their batch is returned. */
	public static void sequences(final String path, final Consumer<Read> consumer){
		Throwable failure=null;
		try(SequenceInput in=new SequenceInput(StreamerFactory.makeStreamer(FileFormat.testInput(path, FileFormat.FASTA, null, true, true), null, true, -1, false, true, 1))){
			in.stream.start();
			for(ListNum<Read> batch=in.stream.nextList(); batch!=null; batch=in.stream.nextList()){
				for(final Read read:batch.list){
					if(failure!=null){continue;}
					try{
						require(read!=null && read.mate==null && read.id!=null && read.bases!=null, "Saved sequence inputs require intact unpaired FASTA records");
						consumer.accept(read);
					}catch(RuntimeException | Error invalid){failure=invalid;}
				}
				in.stream.returnList(batch);
			}
			if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
			if(failure instanceof Error){throw (Error)failure;}
		}
	}
	private void parse(final byte[] row, final int terms){
		p.set(row);
		require(p.terms()==terms, "Saved cache schema changed; expected fields="+terms);
	}
	private static String disposition(final int compatible, final int unique, final int rank){
		return compatible==0 ? "NO_COMPATIBLE_MODEL" : unique==0 ? "AMBIGUOUS_SPANS_ONLY"
			: rank<0 ? "ZERO_NATIVE_INDEX_SUPPORT" : "RANKED";
	}
	static void header(final Input in, final String expected){
		final byte[] row=in.reader.nextLine();
		require(row!=null && Arrays.equals(row, expected.getBytes(StandardCharsets.UTF_8)), "Saved cache header must match the accepted producer schema");
	}
	static String text(final LineParser1 p, final int field){
		p.setBounds(field);
		final String value=new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
		require(value.indexOf('\ufffd')<0, "Invalid UTF-8 cannot change saved input or parent identity");
		return value;
	}
	private static boolean sameFields(final LineParser1 a, final int af, final LineParser1 b, final int bf){
		a.setBounds(af);
		b.setBounds(bf);
		if(a.b()-a.a()!=b.b()-b.a()){return false;}
		for(int i=0; i<a.b()-a.a(); i++){if(a.line()[a.a()+i]!=b.line()[b.a()+i]){return false;}}
		return true;
	}
	private static boolean same(final LineParser1 p, final int field, final byte[] value){
		p.setBounds(field);
		if(p.b()-p.a()!=value.length){return false;}
		for(int i=0; i<value.length; i++){if(p.line()[p.a()+i]!=value[i]){return false;}}
		return true;
	}
	private static boolean bool(final LineParser1 p, final int field){
		require(p.termEquals("true", field) || p.termEquals("false", field), "Saved boolean fields must retain explicit true/false values");
		return p.termEquals("true", field);
	}
	private static int stratum(final LineParser1 p){
		for(int i=0; i<STRATA.length; i++){if(p.termEquals(STRATA[i], 3)){return i;}}
		throw new IllegalArgumentException("Saved selected parents require explicit NOT_SEEN, SEEN or UNKNOWN labels");
	}
	static int minimum(final int a, final int b){return a<0 ? b : b<0 ? a : Math.min(a, b);}
	static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Retained State       ----------------*/
	/*--------------------------------------------------------------*/

	static final class Panel{
		Panel(final String name_, final String prefix_, final String labels_, final int inputCount, final int selected_, final int targets_, final int common_, final int zero_, final int gaps_){
			require(inputCount>0 && selected_>0 && selected_<=inputCount && targets_>=0 && common_>=0 && zero_>=0 && gaps_>=0, "Each saved panel requires complete accepted source counts");
			name=name_;
			nameBytes=name.getBytes(StandardCharsets.UTF_8);
			prefix=prefix_;
			labels=labels_;
			inputs=new Selected[inputCount];
			expectedSelected=selected_;
			expectedTargets=targets_;
			expectedCommon=common_;
			expectedZero=zero_;
			expectedGaps=gaps_;
		}
		final String name, prefix, labels;
		final byte[] nameBytes;
		final Selected[] inputs;
		final int expectedSelected, expectedTargets, expectedCommon, expectedZero, expectedGaps;
		int selected;
		final ObjectMap<String, Parent> parents=new ObjectMap<String, Parent>(String.class, Parent.class);
		final ObjectMap<String, Window> windows=new ObjectMap<String, Window>(String.class, Window.class);
	}
	static final class Selected{
		Selected(final int index_, final Panel panel_, final int slot_, final String id_, final String path_, final int group_, final LineParser1 row){
			index=index_;
			panel=panel_;
			slot=slot_;
			id=id_;
			path=path_;
			group=group_;
			expectedTargets=row.parseInt(5);
			expectedCommon=row.parseInt(6);
			expectedZero=row.parseInt(7);
			expectedGaps=row.parseInt(8);
			expectedNoWindow=row.parseInt(9);
			expectedAmbiguous=row.parseInt(10);
			require(expectedTargets>=0 && expectedCommon>=0 && expectedZero>=0 && expectedGaps>=0
				&& expectedNoWindow>=0 && expectedAmbiguous>=0
				&& (long)expectedCommon+expectedZero+expectedGaps+expectedNoWindow+expectedAmbiguous==expectedTargets,
				"Selected input dispositions must partition all targets, preserving zero-target groups");
			idBytes=id.getBytes(StandardCharsets.UTF_8);
			pathBytes=path.getBytes(StandardCharsets.UTF_8);
		}
		final Panel panel;
		final String id, path;
		final byte[] idBytes, pathBytes;
		final int index, slot, group, expectedTargets, expectedCommon, expectedZero, expectedGaps, expectedNoWindow, expectedAmbiguous;
		int parents, common, zero, gaps, noWindow, ambiguous;
		long expectedWindows, windows;
	}
	static final class Parent{
		Parent(final int index_, final Selected input_, final String id, final String seqid, final int start_, final int stop_, final String strand_, final int stratum_){
			require(start_>0 && stop_>=start_ && (strand_.equals("+") || strand_.equals("-")), "Saved selected parents require valid one-based inclusive bounds and strand");
			index=index_;
			input=input_;
			idBytes=id.getBytes(StandardCharsets.UTF_8);
			seqidBytes=seqid.getBytes(StandardCharsets.UTF_8);
			start=start_;
			stop=stop_;
			strand=strand_.equals("+") ? 0 : 1;
			stratum=stratum_;
			Arrays.fill(controlMin, -1);
			Arrays.fill(observedControls, -1);
		}
		final int index, start, stop, strand, stratum;
		String reason;
		final Selected input;
		final byte[] idBytes, seqidBytes;
		final int[] controlMin=new int[CONTROLS], observedControls=new int[CONTROLS];
		int nativeRank=-1, observedNative=-1, windows, compatible, unique, ambiguous;
		int events, compatibleEvents, uniqueEvents, ambiguousEvents;
		boolean locusRead;
	}
	static final class Window{
		Window(final Selected input_, final long event_, final long ordinal_, final String seqid, final int strand_, final int start_, final int stop_, final int contig, final int length_){
			require(event_>=0 && ordinal_>=0 && (strand_==0 || strand_==1) && start_>=0 && stop_>=start_ && stop_<contig && stop_-start_+1==length_, "Saved windows require intact inclusive bounds, event identity and orientation");
			input=input_;
			event=event_;
			ordinal=ordinal_;
			strand=strand_;
			length=length_;
			start=start_;
			stop=stop_;
			seqidBytes=seqid.getBytes(StandardCharsets.UTF_8);
			key=input.slot+"_"+event;
		}
		final Selected input;
		final long event, ordinal;
		final int strand, length, start, stop;
		final String key;
		final byte[] seqidBytes;
		byte[] bases;
		char[] nativeCounts;
		long[] acceptedModelBits;
		double gc=Double.NaN;
		int nativeModels, acceptedModels;
		int acceptedStart=Integer.MAX_VALUE, acceptedStop=-1;
		final ArrayList<Edge> edges=new ArrayList<Edge>();
	}
	static final class Edge{
		Edge(final Window window_, final Parent parent_, final int models){
			assert(parent_.input==window_.input) : "A cached event cannot acquire a parent from another original input";
			window=window_;
			parent=parent_;
			seenModels=new long[(models+63)/64];
			uniqueModels=new long[seenModels.length];
			Arrays.fill(controls, Integer.MIN_VALUE);
		}
		final Window window;
		final Parent parent;
		final long[] methods=new long[(METHODS+63)/64], seenModels, uniqueModels;
		final int[] controls=new int[CONTROLS*4];
		final int[] nativeTuple=new int[4];
		int methodCount, compatible, unique, nativeRank=-1, links, uniqueLinks;
	}

	final ArrayList<Panel> panels=new ArrayList<Panel>();
	final ArrayList<Selected> selected=new ArrayList<Selected>();
	final ArrayList<Parent> parents=new ArrayList<Parent>();
	final ArrayList<Window> windows=new ArrayList<Window>();
	final int models;
	private final ObjectMap<String, Panel> byName=new ObjectMap<String, Panel>(String.class, Panel.class);
	private final LineParser1 p=new LineParser1('\t');
	private Window cachedWindow;
	private Edge cachedEdge;
	private boolean nativeStarted;
	static final int METHODS=148, CONTROLS=16;
	static final String[] STRATA={"NOT_SEEN", "SEEN", "UNKNOWN"};
	static final int[] CONTROL_METHODS={52,54,55,57,58,60,61,63,73,75,76,78,79,81,82,84};
	private static final int[] CONTROL_INDEX=new int[METHODS];
	static{
		Arrays.fill(CONTROL_INDEX, -1);
		for(int i=0; i<CONTROL_METHODS.length; i++){CONTROL_INDEX[CONTROL_METHODS[i]]=i;}
	}
	static final String INVENTORY_HEADER="panel\tinput_slot\tinput_id\tfasta\tresolved_assembly_group\ttarget_loci\tcommon_loci\tnative_zero_support\tno_compatible_model\tno_window\tambiguous_only";
	static final String RANK_HEADER="input_slot\tinput_id\tevent\tcontig_ordinal\tseqid\tparent\tmethod\twindow_length\tlocus_length\tcompatible_models\tunique_models\trank\ttie_first\ttie_last\tmodel_index\treason";
	static final String LINK_HEADER="input_slot\tinput_id\tevent\tcontig_ordinal\tseqid\tmodel_index\tparent\ttarget\tspan_parent_links\tuniquely_attributed";
	static final String ORACLE_HEADER="input_slot\tinput_id\tevent\tcontig_ordinal\tseqid\tmodel_index\tmodel\tnative_shared_keys\tnative_shortlist_rank\tfully_accepted\taccepted_start1\taccepted_stop1\tpartial\tscore\tscored_parent_links";
}
