package ml;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.PreParser;
import shared.Tools;
import structures.ByteBuilder;
import structures.FloatList;
import structures.IntList;

/**
 * Grows a .bbnet network and optionally deletes low-magnitude active edges.
 * Existing neurons retain their layer position and retained weights exactly.
 *
 * @author Yelan
 */
public final class ModifyNN {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, ModifyNN.class, false);
		new ModifyNN(pp.args, pp.outstream).process();
	}

	ModifyNN(String[] args, PrintStream outstream_){
		outstream=outstream_;
		for(String arg:args){
			final String[] split=arg.split("=", 2);
			final String key=split[0].toLowerCase();
			final String value=split.length>1 ? split[1] : null;
			if(key.equals("in") || key.equals("netin")){in=value;}
			else if(key.equals("out") || key.equals("netout")){out=value;}
			else if(key.equals("dims")){dimsText=value;}
			else if(key.equals("seed")){seed=Long.parseLong(value);}
			else if(key.equals("newweight")){newWeight=Float.parseFloat(value);}
			else if(key.equals("pruneabs")){pruneAbs=Float.parseFloat(value);}
			else if(key.equals("zero2epsilon") || key.equals("zeros2epsilon") || key.equals("zeroes2epsilon") ||
					key.equals("zero2eps")){
				zero2epsilon=Parse.parseBoolean(value);
			}
			else if(key.equals("report")){report=value;}
			else if(key.equals("overwrite") || key.equals("ow")){overwrite=Parse.parseBoolean(value);}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		if(in==null || out==null){throw new IllegalArgumentException(USAGE);}
		Tools.testForDuplicateFiles(true, in, out, report);
		if(!Tools.testOutputFiles(overwrite, false, false, out, report)){
			throw new IllegalArgumentException("Output exists and overwrite=f");
		}
		checkFiniteNonnegative(pruneAbs, "pruneabs");
		if(seed<0){throw new IllegalArgumentException("seed must be nonnegative: "+seed);}
		if(!Float.isFinite(newWeight) || newWeight<MIN_NEW_WEIGHT || newWeight>MAX_NEW_WEIGHT){
			throw new IllegalArgumentException("newweight must be finite and FP16-safe, in ["+
				MIN_NEW_WEIGHT+", "+MAX_NEW_WEIGHT+"]: "+newWeight);
		}
	}

	private void process(){
		final byte[] parentBytes=readFile(in);
		final String parentSha80=sha80(parentBytes);
		final File parentSnapshot=tempInputFile(in, out);
		try{
			writeFile(parentBytes, parentSnapshot);
			final CellNet net=CellNetParser.load(parentSnapshot.getPath(), false);
			validate(net);
			final int[] oldDims=net.dims.clone();
			final int[] newDims=parseDims(dimsText, oldDims);
			final boolean noOp=Arrays.equals(oldDims, newDims) && pruneAbs==0 && !zero2epsilon;
			final Report r=new Report(oldDims, newDims, parentSha80);
			if(noOp){
				checkMatchingTransport(in, out);
				countExistingEdges(net, r);
				publishBytes(parentBytes, out);
				if(report!=null){publishReport(r);}
				outstream.println("ModifyNN copied validated no-op "+in+" -> "+out);
				return;
			}

			final CellNet grown=grow(net, newDims, r);
			grown.setWeightBits(32);
			CellNet.codingA48Out=true;
			publishNet(grown, out);
			if(report!=null){publishReport(r);}
			outstream.println("ModifyNN wrote "+in+" -> "+out+" dims="+Arrays.toString(newDims)+
				" added="+r.added+" deleted="+r.deleted+" converted="+r.converted);
		}finally{
			deleteStagingDirectory(parentSnapshot);
		}
	}

	private CellNet grow(final CellNet old, final int[] newDims, final Report report){
		if(hasPartitionedOutputs(old)){
			throw new IllegalArgumentException("Cannot modify output_partition networks until explicit layout preservation is implemented");
		}
		assert(newDims.length==old.dims.length) : "grow preserves layer count validated by parseDims: new="+
			dimsString(newDims)+" old="+dimsString(old.dims);
		final boolean sourceDense=old.net[1][0].inputs==null;
		final boolean sourceZeroAbsent=sourceDense && old.zeroAbsent();
		final boolean outputDense=sourceDense && pruneAbs==0;
		CellNet.DENSE=outputDense;
		final CellNet grown=new CellNet(newDims, old.seed, old.density, old.density1,
			old.edgeBlockSize, new ArrayList<String>(old.commands));
		grown.cutoff=old.cutoff;
		grown.setZeroAbsent(outputDense && sourceZeroAbsent);
		grown.tags=new LinkedHashMap<String,String>(old.tags);
		archiveParentStats(old, grown);
		grown.tags.put("parent_sha80", report.parentSha80);
		archiveStaleChildTags(grown);
		grown.commands.add("#CL modifynn.sh seed="+seed+" newweight="+newWeight+" pruneabs="+pruneAbs+
			" zero2epsilon="+zero2epsilon+" olddims="+dimsString(old.dims)+" newdims="+dimsString(newDims));
		copyNormalization(old, grown);
		final EdgeRandom random=new EdgeRandom(seed);
		for(int layer=1; layer<newDims.length; layer++){
			final Cell[] oldLayer=old.net[layer], newLayer=grown.net[layer];
			final int oldPrev=old.dims[layer-1], newPrev=newDims[layer-1];
			Function function=null;
			for(int sink=0; sink<newLayer.length; sink++){
				final Cell dest=newLayer[sink];
				final Cell source=sink<oldLayer.length ? oldLayer[sink] : null;
				if(source==null){
					if(function==null){function=unambiguousFunction(oldLayer);}
					dest.function=function;
					dest.bias=0;
				}else{
					dest.function=source.function;
					dest.bias=source.bias;
				}
				if(outputDense){fillDense(source, dest, oldPrev, newPrev, sink, oldLayer.length, random, report, layer, sourceZeroAbsent);}
				else if(sourceDense){fillDenseToSparse(source, dest, oldPrev, newPrev, sink, oldLayer.length, random, report, layer, sourceZeroAbsent);}
				else{fillSparse(source, dest, oldPrev, newPrev, sink, oldLayer.length, random, report, layer);}
			}
		}
		if(!outputDense){CellNet.makeOutputSets(grown.net);}
		grown.makeWeightMatrices();
		countDisconnected(grown, report);
		return grown;
	}

	/** Emits a dense destination row, preserving old columns in place and appending new columns. */
	private void fillDense(Cell source, Cell dest, int oldPrev, int newPrev, int sink, int oldSinks,
			EdgeRandom random, Report report, int layer, boolean sourceZeroAbsent){
		assert(dest.inputs==null) : "fillDense requires the dense destination selected by grow";
		assert(source==null || source.inputs==null) : "fillDense requires a validated dense parent row or a new sink";
		dest.weights=new float[newPrev];
		dest.deltas=new float[newPrev];
		for(int input=0; input<newPrev; input++){
			if(source!=null && input<oldPrev){
				final float w=source.weights[input];
				if(w==0){
					if(sourceZeroAbsent){dest.weights[input]=0;}
					else if(zero2epsilon){dest.weights[input]=newRandomWeight(random); report.converted(layer);}
					else{dest.weights[input]=0; report.retained(layer);}
				}
				else if(Math.abs(w)<pruneAbs){report.deleted(layer); dest.weights[input]=0;}
				else{dest.weights[input]=w; report.retained(layer);}
			}else if(sink>=oldSinks || input>=oldPrev){
				dest.weights[input]=newRandomWeight(random);
				report.added(layer);
			}
		}
	}

	/** Emits a sparse destination row from a sparse parent row, preserving sorted old input indexes. */
	private void fillSparse(Cell source, Cell dest, int oldPrev, int newPrev, int sink, int oldSinks,
			EdgeRandom random, Report report, int layer){
		assert(source==null || source.inputs!=null) : "fillSparse requires a validated sparse parent row or a new sink";
		final IntList inputs=new IntList(newPrev);
		final FloatList weights=new FloatList(newPrev);
		if(source!=null){
			for(int i=0; i<source.inputs.length; i++){
				final int input=source.inputs[i];
				if(input<0 || input>=oldPrev){throw new IllegalArgumentException("Bad sparse input index: "+input);}
				final float w=source.weights[i];
				if(w==0){
					if(zero2epsilon){inputs.add(input); weights.add(newRandomWeight(random)); report.converted(layer);}
					else{inputs.add(input); weights.add(0); report.retained(layer);}
				}
				else if(Math.abs(w)<pruneAbs){report.deleted(layer);}
				else{inputs.add(input); weights.add(w); report.retained(layer);}
			}
		}
		if(sink>=oldSinks){
			for(int input=0; input<newPrev; input++){inputs.add(input); weights.add(newRandomWeight(random)); report.added(layer);}
		}else{
			for(int input=oldPrev; input<newPrev; input++){inputs.add(input); weights.add(newRandomWeight(random)); report.added(layer);}
		}
		dest.inputs=inputs.toArray();
		dest.weights=new float[weights.size()];
		for(int i=0; i<weights.size(); i++){dest.weights[i]=weights.get(i);}
		dest.deltas=new float[dest.weights.length];
	}

	/** Emits an explicit sparse row from a dense parent row so pruned dense edges vanish. */
	private void fillDenseToSparse(Cell source, Cell dest, int oldPrev, int newPrev, int sink, int oldSinks,
			EdgeRandom random, Report report, int layer, boolean sourceZeroAbsent){
		assert(source==null || source.inputs==null) : "fillDenseToSparse indexes a dense parent row before deleting edges";
		final IntList inputs=new IntList(newPrev);
		final FloatList weights=new FloatList(newPrev);
		if(source!=null){
			for(int input=0; input<oldPrev; input++){
				final float w=source.weights[input];
				if(w==0){
					if(!sourceZeroAbsent){
						if(zero2epsilon){inputs.add(input); weights.add(newRandomWeight(random)); report.converted(layer);}
						else{inputs.add(input); weights.add(0); report.retained(layer);}
					}
				}
				else if(Math.abs(w)<pruneAbs){report.deleted(layer);}
				else{inputs.add(input); weights.add(w); report.retained(layer);}
			}
		}
		if(sink>=oldSinks){
			for(int input=0; input<newPrev; input++){inputs.add(input); weights.add(newRandomWeight(random)); report.added(layer);}
		}else{
			for(int input=oldPrev; input<newPrev; input++){inputs.add(input); weights.add(newRandomWeight(random)); report.added(layer);}
		}
		dest.inputs=inputs.toArray();
		dest.weights=new float[weights.size()];
		for(int i=0; i<weights.size(); i++){dest.weights[i]=weights.get(i);}
		dest.deltas=new float[dest.weights.length];
	}

	private float newRandomWeight(EdgeRandom random){
		return random.nextWeight(newWeight);
	}

	private static void archiveParentStats(final CellNet old, final CellNet grown){
		grown.tags.put("parent_error_rate", Float.toString(old.errorRate));
		grown.tags.put("parent_weighted_error_rate", Float.toString(old.weightedErrorRate));
		grown.tags.put("parent_fp_rate", Float.toString(old.fpRate));
		grown.tags.put("parent_fn_rate", Float.toString(old.fnRate));
		grown.tags.put("parent_epoch", Integer.toString(old.epoch));
		grown.tags.put("parent_epochs_trained", Long.toString(old.epochsTrained));
		grown.tags.put("parent_samples_trained", Long.toString(old.samplesTrained));
		if(old.lastStats!=null){grown.tags.put("parent_last_stats", old.lastStats);}
	}

	private static void archiveStaleChildTags(final CellNet grown){
		final LinkedHashMap<String,String> archive=new LinkedHashMap<String,String>();
		for(Map.Entry<String,String> entry : grown.tags.entrySet()){
			final String key=entry.getKey().toLowerCase();
			if(staleChildTag(key)){archive.put("parent_"+key, entry.getValue());}
		}
		if(archive.isEmpty()){return;}
		final ArrayList<String> remove=new ArrayList<String>();
		for(String key : grown.tags.keySet()){
			if(staleChildTag(key.toLowerCase())){remove.add(key);}
		}
		for(String key : remove){grown.tags.remove(key);}
		grown.tags.putAll(archive);
	}

	private static boolean staleChildTag(final String key){
		if(key.startsWith("magqc_")){return true;}
		for(String stale : STALE_CHILD_TAGS){
			if(key.equals(stale)){return true;}
		}
		return false;
	}

	/** Copies explicit input standardization and appends identity-normalized columns. */
	private void copyNormalization(final CellNet old, final CellNet grown){
		final float[] oldMean=old.inputMeanCopy(), oldInverseStd=old.inputInverseStdCopy();
		if(oldMean==null){return;}
		assert(oldMean.length==old.dims[0] && oldInverseStd.length==old.dims[0]) :
			"Normalization copies must cover the input width validated by CellNet: inputs="+old.dims[0]+
			" means="+oldMean.length+" inverseStd="+oldInverseStd.length;
		assert(grown.dims[0]>=old.dims[0]) : "copyNormalization extends inputs without truncating parent statistics: new="+
			grown.dims[0]+" old="+old.dims[0];
		final float[] mean=new float[grown.dims[0]], inverseStd=new float[grown.dims[0]];
		System.arraycopy(oldMean, 0, mean, 0, oldMean.length);
		System.arraycopy(oldInverseStd, 0, inverseStd, 0, oldInverseStd.length);
		Arrays.fill(inverseStd, oldInverseStd.length, inverseStd.length, 1f);
		grown.setInputNormalization(mean, inverseStd);
	}

	/** Returns a single old activation for new neurons, rejecting mixed activations in grown layers. */
	private static Function unambiguousFunction(final Cell[] layer){
		Function f=layer[0].function;
		for(Cell cell:layer){
			if(cell.function!=f){
				throw new IllegalArgumentException("New neurons require an unambiguous old activation per layer");
			}
		}
		return f;
	}

	private static int[] parseDims(String text, int[] oldDims){
		if(text==null || text.equalsIgnoreCase("null")){return oldDims.clone();}
		String[] split=text.split(",");
		if(split.length!=oldDims.length){
			throw new IllegalArgumentException("dims arity "+split.length+" != old layer count "+oldDims.length);
		}
		int[] dims=new int[split.length];
		for(int i=0; i<dims.length; i++){
			dims[i]=Integer.parseInt(split[i]);
			if(dims[i]<oldDims[i]){
				throw new IllegalArgumentException("dims may not shrink layer "+i+": "+dims[i]+" < "+oldDims[i]);
			}
		}
		return dims;
	}

	private static void checkFiniteNonnegative(float value, String name){
		if(!Float.isFinite(value) || value<0){throw new IllegalArgumentException(name+" must be finite and nonnegative");}
	}

	private static String dimsString(int[] dims){
		ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<dims.length; i++){
			if(i>0){bb.append(',');}
			bb.append(dims[i]);
		}
		return bb.toString();
	}

	private static boolean hasPartitionedOutputs(final CellNet net){
		return net.getTag("output_partition")!=null;
	}

	/** Validates finite weights and row topology before a no-op copy or in-place growth. */
	private static void validate(final CellNet net){
		final boolean dense=net.net[1][0].inputs==null;
		for(int layer=1; layer<net.net.length; layer++){
			final int prev=net.dims[layer-1];
			for(Cell cell:net.net[layer]){
				if(cell==null || cell.function==null){throw new IllegalArgumentException("Missing cell at layer "+layer);}
				if(!Float.isFinite(cell.bias())){throw new IllegalArgumentException("Nonfinite bias at cell "+cell.id());}
				if(cell.weights==null){throw new IllegalArgumentException("Missing weights at cell "+cell.id());}
				if(cell.deltas!=null && cell.deltas.length!=cell.weights.length){
					throw new IllegalArgumentException("Delta/weight length mismatch at cell "+cell.id());
				}
				if(dense){
					if(cell.inputs!=null || cell.weights.length!=prev){
						throw new IllegalArgumentException("Dense topology mismatch at cell "+cell.id());
					}
				}else{
					if(cell.inputs==null || cell.inputs.length!=cell.weights.length){
						throw new IllegalArgumentException("Sparse topology mismatch at cell "+cell.id());
					}
					int last=-1;
					for(int input : cell.inputs){
						if(input<=last || input<0 || input>=prev){
							throw new IllegalArgumentException("Sparse input index mismatch at cell "+cell.id());
						}
						last=input;
					}
				}
				for(float weight : cell.weights){
					if(!Float.isFinite(weight)){throw new IllegalArgumentException("Nonfinite weight at cell "+cell.id());}
				}
			}
		}
	}

	/** Counts all active parent edges for a validated no-op copy report. */
	private static void countExistingEdges(final CellNet net, final Report report){
		final boolean dense=net.net[1][0].inputs==null;
		final boolean zeroAbsent=dense && net.zeroAbsent();
		for(int layer=1; layer<net.net.length; layer++){
			for(Cell cell:net.net[layer]){
				final int count=countLiveIncoming(cell, dense, zeroAbsent);
				for(int i=0; i<count; i++){report.retained(layer);}
			}
		}
		countDisconnected(net, report);
	}

	private static void countDisconnected(final CellNet net, final Report report){
		final boolean dense=net.net[1][0].inputs==null;
		final boolean zeroAbsent=dense && net.zeroAbsent();
		for(int layer=1; layer<net.net.length; layer++){
			for(Cell cell:net.net[layer]){
				if(countLiveIncoming(cell, dense, zeroAbsent)==0){report.disconnected[layer-1]++;}
			}
		}
	}

	/** Counts structurally live incoming weights under the parsed dense zero-mask contract. */
	private static int countLiveIncoming(final Cell cell, final boolean dense, final boolean zeroAbsent){
		int count=0;
		final int limit=dense ? cell.weights.length : cell.inputs.length;
		for(int i=0; i<limit; i++){
			if(!dense || !zeroAbsent || cell.weights[i]!=0){count++;}
		}
		return count;
	}

	private void publishNet(CellNet net, String path){
		final File temp=tempFile(path);
		boolean moved=false;
		try{
			writeNet(net, temp.getPath());
			moveTemp(temp, path);
			moved=true;
		}finally{
			if(!moved){deleteStagingDirectory(temp);}
		}
	}

	/** Publishes fully written bytes by rename from a private directory beside the destination. */
	private void publishBytes(byte[] bytes, String dest){
		final File temp=tempFile(dest);
		boolean moved=false;
		try{
			writeFile(bytes, temp);
			moveTemp(temp, dest);
			moved=true;
		}finally{
			if(!moved){deleteStagingDirectory(temp);}
		}
	}

	/** Publishes the TSV report with the same checked writer/rename path as .bbnet output. */
	private void publishReport(Report r){
		final File temp=tempFile(report);
		boolean moved=false;
		try{
			writeReport(r, temp.getPath());
			moveTemp(temp, report);
			moved=true;
		}finally{
			if(!moved){deleteStagingDirectory(temp);}
		}
	}

	private static File tempFile(String path){
		return tempFile(path, new File(path).getAbsoluteFile().getName());
	}

	private static File tempInputFile(String inputPath, String outputPath){
		return tempFile(outputPath, new File(inputPath).getAbsoluteFile().getName());
	}

	/** Allocates a private staging path under the final output directory. */
	private static File tempFile(String path, String stagedName){
		final File target=new File(path).getAbsoluteFile();
		final File parent=target.getParentFile();
		if(parent==null || !parent.isDirectory()){
			throw new IllegalArgumentException("Output parent is not a directory: "+parent);
		}
		try{
			final Path dir=Files.createTempDirectory(parent.toPath(), ".modifynn.");
			return dir.resolve(stagedName).toFile();
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}

	/** Renames a finished staged file into place, honoring overwrite only at publication. */
	private void moveTemp(final File temp, final String path){
		final File target=new File(path);
		try{
			if(overwrite){Files.move(temp.toPath(), target.toPath(), StandardCopyOption.REPLACE_EXISTING);}
			else{Files.move(temp.toPath(), target.toPath());}
		}catch(IOException e){
			throw new RuntimeException(e);
		}finally{
			deleteStagingDirectory(temp);
		}
	}

	/** Deletes only the private directory family created by tempFile. */
	private static void deleteStagingDirectory(final File temp){
		final File dir=temp.getParentFile();
		if(dir==null){return;}
		final String name=dir.getName();
		if(!name.startsWith(".modifynn.")){return;}
		final File[] children=dir.listFiles();
		if(children!=null){for(File child : children){child.delete();}}
		dir.delete();
	}

	/** Reads a stable raw byte image for exact no-op copies and parent provenance. */
	private static byte[] readFile(final String path){
		try{
			final FileInputStream in=new FileInputStream(path);
			try{
				final FileChannel channel=in.getChannel();
				final long size=channel.size();
				if(size>Integer.MAX_VALUE){throw new IllegalArgumentException("Input is too large: "+path);}
				final byte[] array=new byte[(int)size];
				int off=0;
				for(int len; off<array.length && (len=in.read(array, off, array.length-off))>=0; ){
					off+=len;
				}
				if(off<array.length){throw new IOException("Short read from "+path);}
				if(in.read()>=0){throw new IOException("Input grew while reading "+path);}
				return array;
			}finally{in.close();}
		}catch(IOException e){throw new RuntimeException(e);}
	}

	/** Rejects no-op copies that would change compression by raw byte transport. */
	private static void checkMatchingTransport(String a, String b){
		final String at=ReadWrite.compressionType(a), bt=ReadWrite.compressionType(b);
		if(at==null ? bt!=null : !at.equals(bt)){
			throw new IllegalArgumentException("No-op byte-copy requires matching input/output compression suffixes: "+a+", "+b);
		}
	}

	/** Writes bytes synchronously into a private staging file. */
	private static void writeFile(final byte[] bytes, final File file){
		try{
			final FileOutputStream out=new FileOutputStream(file);
			try{out.write(bytes);}
			finally{out.close();}
		}catch(IOException e){throw new RuntimeException(e);}
	}

	/** Returns the last 80 bits of the SHA-256 digest as 20 lowercase hex characters. */
	private static String sha80(final byte[] bytes){
		try{
			final byte[] digest=MessageDigest.getInstance("SHA-256").digest(bytes);
			final char[] out=new char[20];
			for(int i=0; i<10; i++){
				final int x=digest[digest.length-10+i]&0xff;
				out[2*i]=HEX[x>>>4]; out[2*i+1]=HEX[x&15];
			}
			return new String(out);
		}catch(Exception e){throw new RuntimeException(e);}
	}

	private static void writeNet(CellNet net, String path){
		FileFormat ff=FileFormat.testOutput(path, FileFormat.BBNET, null, true, true, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.println(net.toBytes());
		if(bsw.poisonAndWait()){throw new IllegalStateException("Could not write "+path);}
	}

	/** Writes a checked TSV report; ByteStreamWriter preserves gzip suffix handling and reports close errors. */
	private void writeReport(Report r, String path){
		final ByteBuilder bb=new ByteBuilder();
		bb.append("key\tvalue\n");
		bb.append("old_dims\t").append(dimsString(r.oldDims)).nl();
		bb.append("new_dims\t").append(dimsString(r.newDims)).nl();
		bb.append("seed\t").append(seed).nl();
		bb.append("parent_sha80\t").append(r.parentSha80).nl();
		bb.append("newweight\t").append(Float.toString(newWeight)).nl();
		bb.append("pruneabs\t").append(Float.toString(pruneAbs)).nl();
		bb.append("zero2epsilon\t").append(zero2epsilon).nl();
		bb.append("index_preservation\tappend_only").nl();
		bb.append("retained_edges\t").append(r.retained).nl();
		bb.append("added_edges\t").append(r.added).nl();
		bb.append("deleted_edges\t").append(r.deleted).nl();
		bb.append("converted_edges\t").append(r.converted).nl();
		for(int i=0; i<r.retainedByLayer.length; i++){
			final int layer=i+1;
			bb.append("layer_").append(layer).append("_retained_edges\t").append(r.retainedByLayer[i]).nl();
			bb.append("layer_").append(layer).append("_added_edges\t").append(r.addedByLayer[i]).nl();
			bb.append("layer_").append(layer).append("_deleted_edges\t").append(r.deletedByLayer[i]).nl();
			bb.append("layer_").append(layer).append("_converted_edges\t").append(r.convertedByLayer[i]).nl();
			bb.append("layer_").append(layer).append("_disconnected_neurons\t").append(r.disconnected[i]).nl();
		}
		FileFormat ff=FileFormat.testOutput(path, FileFormat.TXT, null, true, true, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.print(bb);
		if(bsw.poisonAndWait()){throw new IllegalStateException("Could not write "+path);}
	}

	/** Per-run edge accounting for the optional TSV report. */
	private static final class Report{
		Report(int[] oldDims_, int[] newDims_, String parentSha80_){
			oldDims=oldDims_; newDims=newDims_; parentSha80=parentSha80_;
			retainedByLayer=new long[oldDims.length-1];
			addedByLayer=new long[oldDims.length-1];
			deletedByLayer=new long[oldDims.length-1];
			convertedByLayer=new long[oldDims.length-1];
			disconnected=new long[oldDims.length-1];
		}
		void retained(int layer){retained++; retainedByLayer[layer-1]++;}
		void added(int layer){added++; addedByLayer[layer-1]++;}
		void deleted(int layer){deleted++; deletedByLayer[layer-1]++;}
		void converted(int layer){converted++; convertedByLayer[layer-1]++;}
		final int[] oldDims, newDims;
		final String parentSha80;
		final long[] retainedByLayer, addedByLayer, deletedByLayer, convertedByLayer, disconnected;
		long retained, added, deleted, converted;
	}

	/** Deterministic SplitMix64 stream for appended edge weights. */
	private static final class EdgeRandom{
		EdgeRandom(long seed_){state=seed_;}

		float nextWeight(final float epsilon){
			final long word=Tools.splitMix64(state+=GAMMA);
			final float u=(word&0xFFFFFFL)*0x1p-24f;
			final float magnitude=epsilon*(0.5f+0.5f*u);
			return word<0 ? -magnitude : magnitude;
		}

		private long state;
	}

	public static final String USAGE="modifynn.sh in=old.bbnet out=grown.bbnet "+
		"dims=N,H1,H2,O seed=1 newweight=1e-3 pruneabs=0 zero2epsilon=f report=changes.tsv";

	private static final long GAMMA=0x9E3779B97F4A7C15L;
	private static final float MIN_NEW_WEIGHT=0x1p-13f;
	private static final float MAX_NEW_WEIGHT=65504f;
	private static final char[] HEX="0123456789abcdef".toCharArray();
	private static final String[] STALE_CHILD_TAGS={
		"affine_prefix_layers", "affine_prefix_semantics", "best_val_mse", "checkpoint_sha256",
		"checkpoint_sha80", "exported_by", "gpu_trainer", "mean_sd_folded", "normalization",
		"output_contract", "output_partition", "repr_note", "training_density", "training_dims"};

	private final PrintStream outstream;
	private String in, out, dimsText, report;
	private boolean overwrite=false, zero2epsilon=false;
	private long seed=1;
	private float newWeight=1e-3f, pruneAbs=0;
}
