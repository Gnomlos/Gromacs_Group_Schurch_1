### Rationale: Hydrated Gd³⁺ and the Role of Fpocket in Binding Site Identification

The interaction of gadolinium ions (Gd³⁺) with proteins critically depends on the shape and size of available protein cavities. In physiological conditions, Gd³⁺ exists not as a bare ion but as a hydrated complex, typically coordinating **8–9 water molecules** in its primary hydration shell :contentReference[oaicite:7]{index=7}. The average Gd–O distance is about **2.37 Å**, implying a hydrated ion diameter of approximately **4.5–5.0 Å** :contentReference[oaicite:8]{index=8}.

<p align="center">
  <!-- Replace these paths with your own images in docs/images -->
  <img src="docs/images/gd3_hydrate_schematic.png" alt="Hydrated Gd³⁺ ion" width="300"/>
</p>

#### Why this matters for using Fpocket

- **Size filtering**: Fpocket's detection of cavities helps pinpoint sites whose volume and geometry can accommodate the hydrated Gd³⁺ (≈ 5 Å diameter), effectively filtering out pockets that are too small.
- **Ranking potential sites**: The tool ranks cavities by volume and physico-chemical descriptors, highlighting those most promising for Gd³⁺ accommodation.
- **Efficiency**: Fpocket serves as a **pre-screening** step in the workflow, enabling focused downstream analyses (e.g., docking or molecular dynamics) on only the most viable binding sites.

In essence, fpocket functions as an **efficient structural filter**, ensuring that only pockets large enough for the hydrated Gd³⁺ complex enter into more computationally intensive investigations.

---

Would you like me to assist with designing a custom schematic (e.g., producing a simple vector diagram of the Gd³⁺ hydration shell with size annotation), or integrating the actual hydration images into your project asset folder?
::contentReference[oaicite:9]{index=9}
