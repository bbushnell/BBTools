package prok;

import java.io.File;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;

import consensus.BaseGraph;

/** Focused text/compressed-text/serialized HBM dispatch fixture. */
public class TrnaConsensusBuilderModelLoadTest {

	public static void main(String[] args) throws Exception{
		final File dir=Files.createTempDirectory("hbm_load_").toFile();
		try{
			final ArrayList<BaseGraph> source=new ArrayList<BaseGraph>();
			source.add(new BaseGraph("universal",bytes("ACGTTGCA"),null,0,0));
			source.add(new BaseGraph("bacteria",bytes("GGGCCCAA"),null,1,0));
			final File hbmGz=new File(dir,"models.hbm.gz"), txtGz=new File(dir,"models.txt.gz"), serialized=new File(dir,"models.ser");
			TrnaConsensusBuilder.writeTextModels(source,hbmGz.getAbsolutePath());
			TrnaConsensusBuilder.writeTextModels(source,txtGz.getAbsolutePath());
			checkModels(TrnaConsensusBuilder.loadModels(hbmGz.getAbsolutePath()),".hbm.gz");
			checkModels(TrnaConsensusBuilder.loadModels(txtGz.getAbsolutePath()),".txt.gz");
			try(ObjectOutputStream out=new ObjectOutputStream(new FileOutputStream(serialized))){out.writeObject(source);}
			checkModels(TrnaConsensusBuilder.loadModels(serialized.getAbsolutePath()),"serialized fallback");
			System.out.println("PASS TrnaConsensusBuilderModelLoadTest");
		}finally{
			final File[] files=dir.listFiles(); if(files!=null){for(File file:files){file.delete();}} dir.delete();
		}
	}

	private static void checkModels(BaseGraph[] models,String label){
		check(models!=null && models.length==2,label+" model count");
		check("universal".equals(models[0].name) && "bacteria".equals(models[1].name),label+" names");
		check(Arrays.equals(bytes("ACGTTGCA"),models[0].original) && Arrays.equals(bytes("GGGCCCAA"),models[1].original),label+" sequences");
	}

	private static byte[] bytes(String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void check(boolean ok,String message){if(!ok){throw new AssertionError(message);}}
}
