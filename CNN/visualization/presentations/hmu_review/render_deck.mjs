import fs from 'node:fs/promises';
import { FileBlob, PresentationFile } from '@oai/artifact-tool';

const root = '/Users/fabioali/SND@LHC/nue_analysis/CNN';
const deck = await PresentationFile.importPptx(
  await FileBlob.load(`${root}/output/hmu_review_en/CNN_SND_Hmu_Review_EN_v4.pptx`),
);
const output = `${root}/visualization/presentations/hmu_review/rendered`;
await fs.mkdir(output, { recursive: true });
for (let i = 0; i < deck.slides.items.length; i += 1) {
  const image = await deck.export({ slide: deck.slides.items[i], format: 'png', scale: 1 });
  await fs.writeFile(
    `${output}/slide-${String(i + 1).padStart(2, '0')}.png`,
    new Uint8Array(await image.arrayBuffer()),
  );
}
console.log(`Rendered ${deck.slides.items.length} slides`);
