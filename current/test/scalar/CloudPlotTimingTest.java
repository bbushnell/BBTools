package test.scalar;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import scalar.CloudPlot;

/** Deterministic real-rendering workload for off-axis color-range performance. */
public class CloudPlotTimingTest{
	/** Generates a modern table, then repeatedly renders it through the production entrypoint. */
	public static void main(String[] args) throws Exception{
		final Path root=Paths.get(args[0]);
		final int rows=Integer.parseInt(args[1]), rounds=Integer.parseInt(args[2]);
		final boolean constant=args.length>3 && args[3].equals("constant");
		if(rows<1 || rounds<1){throw new IllegalArgumentException("rows and rounds must be positive");}
		Files.createDirectories(root);
		final StringBuilder text=new StringBuilder("#Name\tLength\tGC\tHH\tCAGA\tDepth\tStart\tTaxID\tTaxID2\n");
		for(int i=0; i<rows; i++){
			text.append("ctg").append(i).append('\t').append(1000+i%901).append('\t');
			text.append((i%97)/96f).append('\t').append((i%89)/88f).append('\t');
			text.append((i%83)/82f).append('\t').append(1+(i*37)%503).append("\t0\t0\t0\n");
		}
		final Path input=root.resolve("points.tsv"), output=root.resolve("plot.png");
		Files.write(input, text.toString().getBytes(StandardCharsets.UTF_8));
		final long start=System.nanoTime();
		for(int i=0; i<rounds; i++){
			CloudPlot.main(new String[]{"in="+input, "out="+output, "order=gc,hh,caga"+(constant ? ",depth" : ""), "colorby=depth",
					"decorrelate=f", "legend=f", "tree=f", "t=1", "overwrite=t"});
		}
		if(Files.size(output)==0){throw new AssertionError("No rendered PNG");}
		if(constant){
			final java.awt.image.BufferedImage img=javax.imageio.ImageIO.read(output.toFile());
			int colored=0;
			for(int y=0; y<img.getHeight(); y++){
				for(int x=0; x<img.getWidth(); x++){
					final int rgb=img.getRGB(x, y), r=(rgb>>>16)&255, g=(rgb>>>8)&255, b=rgb&255;
					if(r!=g || g!=b){colored++;}
				}
			}
			if(colored==0){throw new AssertionError("Constant-valued point disappeared from the plot");}
		}
		System.out.println("rows="+rows+" rounds="+rounds+" seconds="+((System.nanoTime()-start)*1e-9));
	}
}
