package ml;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;

/** Analytic ownership and topology fixtures for partition-aware native growth.
 * @author Nilou
 */
public final class ModifyNNPartitionTest {

	public static void main(String[] args) throws Exception{
		if(args.length>1){throw new IllegalArgumentException("Expected optional fresh fixture directory");}
		final File root=args.length==0 ? Files.createTempDirectory("modifynn_partition_").toFile() : Files.createDirectory(new File(args[0]).toPath()).toFile();
		try{
			for(boolean dense : new boolean[]{false, true}){checkGrowth(root, dense);}
			checkRejected(root);
		}finally{if(args.length==0){delete(root);}}
		System.out.println("ModifyNNPartitionTest PASS");
	}

	private static void checkGrowth(File root, boolean dense) throws Exception{
		final File in=new File(root, "parent_"+dense+".bbnet"), out=new File(root, "child_"+dense+".bbnet");
		final CellNet parent=parent(dense);
		Files.write(in.toPath(), parent.toBytes().toBytes());
		ModifyNN.main(new String[]{"in="+in, "out="+out, "dims=5,14,4", "privateperhead=3", "seed=1"});
		final CellNet child=CellNetParser.load(out.toString(), false);
		check(Arrays.equals(child.dims, new int[]{5, 14, 4}), "Child dimensions");
		check(child.net[1][0].inputs!=null, "Partition child must retain explicit topology");
		check(Float.floatToRawIntBits(weight(child.net[1][0], 3))==-1174124527, "First appended edge uses shared SplitMix64 contract");
		for(int h=0; h<2; h++){
			check(Arrays.equals(child.net[2][h].inputs, parent.net[2][h].inputs==null ?
				indexes(h) : parent.net[2][h].inputs), "Old head ownership must not acquire new private nodes");
			for(int i : indexes(h)){check(weight(child.net[2][h], i)==weight(parent.net[2][h], i), "Retained final weight changed");}
			check(child.net[2][h].bias==parent.net[2][h].bias, "Old head bias changed");
		}
		check(Arrays.equals(child.net[2][2].inputs, new int[]{0, 1, 2, 3, 8, 9, 10}), "First appended head ownership");
		check(Arrays.equals(child.net[2][3].inputs, new int[]{0, 1, 2, 3, 11, 12, 13}), "Second appended head ownership");
		for(int layer=1; layer<parent.net.length; layer++){
			for(int h=0; h<parent.net[layer].length; h++){
				for(int i=0; i<parent.dims[layer-1]; i++){
					check(weight(child.net[layer][h], i)==weight(parent.net[layer][h], i), "Retained weight changed");
				}
			}
		}
		check(Arrays.equals(child.inputMeanCopy(), new float[]{1, 2, 3, 0, 0}), "Mean prefix preservation");
		check(Arrays.equals(child.inputInverseStdCopy(), new float[]{2, 3, 4, 1, 1}), "Normalization prefix preservation");
		check(child.getTag("output_partition").contains("\"shared_inputs\":[0,1,2,3]"), "Shared ownership metadata");
		check(child.getTag("parent_output_partition").equals(parent.getTag("output_partition")), "Parent metadata archived");
		final File noop=new File(root, "noop_"+dense+".bbnet");
		ModifyNN.main(new String[]{"in="+out, "out="+noop});
		check(Arrays.equals(Files.readAllBytes(out.toPath()), Files.readAllBytes(noop.toPath())), "Partition no-op is byte-identical");
		final File twice=new File(root, "twice_"+dense+".bbnet");
		ModifyNN.main(new String[]{"in="+out, "out="+twice, "dims=5,15,5", "privateperhead=1"});
		final CellNet second=CellNetParser.load(twice.toString(), false);
		for(int h=0; h<4; h++){check(Arrays.equals(second.net[2][h].inputs, child.net[2][h].inputs), "Repeat growth changed old membership");}
		check(Arrays.equals(second.net[2][4].inputs, new int[]{0, 1, 2, 3, 14}), "Repeat growth private ownership");
	}

	private static void checkRejected(File root) throws Exception{
		final File in=new File(root, "reject.bbnet");
		Files.write(in.toPath(), parent(false).toBytes().toBytes());
		reject(in, new File(root, "bad_dims.bbnet"), "dims=5,15,4", "privateperhead=3");
		reject(in, new File(root, "bad_private.bbnet"), "privateperhead=-1");
		final CellNet bad=parent(false);
		// Deliberately corrupt the serialized input index; do not rebuild an already initialized network.
		bad.net[2][0].inputs[5]=6;
		final File forbidden=new File(root, "forbidden.bbnet");
		Files.write(forbidden.toPath(), bad.toBytes().toBytes());
		reject(forbidden, new File(root, "bad_edge.bbnet"), "dims=4,8,2");
	}

	private static void reject(File in, File out, String... options){
		final String[] args=new String[2+options.length];
		args[0]="in="+in; args[1]="out="+out;
		System.arraycopy(options, 0, args, 2, options.length);
		boolean failed=false;
		try{ModifyNN.main(args);}catch(IllegalArgumentException e){failed=true;}
		check(failed && !out.exists(), "Invalid partition request must fail before publication");
	}

	private static CellNet parent(boolean dense){
		CellNet.DENSE=dense;
		final CellNet net=new CellNet(new int[]{3, 8, 2}, 1, 1f, 0f, 1, new ArrayList<String>());
		for(int layer=1; layer<net.net.length; layer++){
			for(int h=0; h<net.net[layer].length; h++){
				final Cell cell=net.net[layer][h];
				final int[] inputs=layer==1 ? new int[]{0, 1, 2} : indexes(h);
				cell.function=Function.getFunction(Function.LINEAR); cell.bias=h/100f;
				cell.inputs=dense ? null : inputs;
				cell.weights=new float[dense ? net.dims[layer-1] : inputs.length];
				for(int j=0; j<inputs.length; j++){cell.weights[dense ? inputs[j] : j]=(h+1)*0.01f+(inputs[j]+1)*0.001f;}
				cell.deltas=new float[cell.weights.length];
			}
		}
		if(!dense){net.net[1][0].weights[0]=0; CellNet.makeOutputSets(net.net);}
		net.setZeroAbsent(dense);
		net.makeWeightMatrices();
		net.setInputNormalization(new float[]{1, 2, 3}, new float[]{2, 3, 4});
		net.setTag("output_partition", "{\"scheme\":\"half_shared_v1\",\"hidden_width\":8,\"outputs\":2,\"shared\":4,\"private_counts\":[2,2],\"head_order\":[0,1]}");
		return net;
	}

	private static int[] indexes(int head){return new int[]{0, 1, 2, 3, 4+2*head, 5+2*head};}

	private static float weight(Cell cell, int input){
		if(cell.inputs==null){return cell.weights[input];}
		final int index=Arrays.binarySearch(cell.inputs, input);
		return index<0 ? 0 : cell.weights[index];
	}

	private static void check(boolean condition, String message){if(!condition){throw new AssertionError(message);}}

	private static void delete(File file){
		if(file.isDirectory()){for(File child : file.listFiles()){delete(child);}}
		if(!file.delete()){throw new IllegalStateException("Could not remove fixture file: "+file);}
	}
}
