package org.panda.resource.tcga;

import org.panda.resource.tcga.MutSigReader;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Ozgun Babur
 */
public class PanMutSigReader
{
	static String dir = "/home/exacloud/lustre1/users/peto/GDAN/";

	public static Map<String, Double> getBestPValues()
	{
		Map<String, Double> best = new HashMap<>();

		for (File d : new File(dir).listFiles())
		{
			if (new File(d.getPath() + "/scores-mutsig.txt").exists())
			{
				Map<String, Double> pvals = MutSigReader.readPValues(d.getPath());

				for (String gene : pvals.keySet())
				{
					if (!best.containsKey(gene) || best.get(gene) > pvals.get(gene))
					{
						best.put(gene, pvals.get(gene));
					}
				}
			}
		}
		return best;
	}

	public static void main(String[] args) throws IOException
	{
		// Write pancan ranking based on mutsig

		Map<String, Double> pvals = getBestPValues();
		BufferedWriter writer = Files.newBufferedWriter(Paths.get("/home/exacloud/lustre1/users/peto/GDAN/RankedGenes.txt"));
		writer.write("Gene\tMutSig");
		pvals.keySet().stream().sorted((g1, g2) -> pvals.get(g1).compareTo(pvals.get(g2))).forEach(g ->
			FileUtil.lnwrite(g + "\t" + pvals.get(g), writer));
		writer.close();
	}
}
