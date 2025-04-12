// Model.cpp
#include "Model.h"
#include <iostream>
#include <glm/gtc/type_ptr.hpp>
#include "stb_image.h"
#include <functional>

Model::Model() :
    currentAnimation(0),
    animationTime(0.0f),
    position(0.0f),
    rotationAngle(0.0f),
    rotationAxis(0.0f, 1.0f, 0.0f),
    scale(1.0f, 1.0f, 1.0f),
    scene(nullptr) {
}

Model::~Model() {
    // Clean up OpenGL resources
    for (auto& mesh : meshes) {
        glDeleteVertexArrays(1, &mesh.VAO);
        glDeleteBuffers(1, &mesh.VBO);
        glDeleteBuffers(1, &mesh.EBO);
    }

    // Clean up textures
    for (auto& mesh : meshes) {
        for (auto& texture : mesh.textures) {
            glDeleteTextures(1, &texture.id);
        }
    }
}

bool Model::loadModel(const std::string& path) {
    // Store the directory for loading associated textures
    directory = path.substr(0, path.find_last_of("/\\"));

    // Import using Assimp
    unsigned int importFlags =
        aiProcess_Triangulate |            // Triangulate polygons (if any)
        aiProcess_GenSmoothNormals |       // Generate normals if not present
        aiProcess_FlipUVs |                // Flip texture coords
        aiProcess_CalcTangentSpace |       // Calculate tangents for normal mapping
        aiProcess_JoinIdenticalVertices |  // Optimize mesh
        aiProcess_LimitBoneWeights |       // Limit bone weights to 4 per vertex
        aiProcess_ImproveCacheLocality;    // Improve cache locality

    scene = importer.ReadFile(path, importFlags);

    // Check for errors
    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
        std::cerr << "ERROR::ASSIMP::" << importer.GetErrorString() << std::endl;
        return false;
    }

    // Process nodes recursively
    processNode(scene->mRootNode, scene);

    // Load animations if available
    loadAnimations(scene);

    // Pre-allocate bone transforms array
    boneTransforms.resize(bones.size(), glm::mat4(1.0f));

    std::cout << "Loaded model: " << path << std::endl;
    std::cout << "  Meshes: " << meshes.size() << std::endl;
    std::cout << "  Bones: " << bones.size() << std::endl;
    std::cout << "  Animations: " << animations.size() << std::endl;

    return true;
}

void Model::processNode(aiNode* node, const aiScene* scene) {
    // Process all meshes in this node
    for (unsigned int i = 0; i < node->mNumMeshes; i++) {
        aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
        meshes.push_back(processMesh(mesh, scene));
    }

    // Process all child nodes recursively
    for (unsigned int i = 0; i < node->mNumChildren; i++) {
        processNode(node->mChildren[i], scene);
    }
}

Model::Mesh Model::processMesh(aiMesh* mesh, const aiScene* scene) {
    Mesh result;

    // Process vertices
    for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
        Vertex vertex;

        // Position
        vertex.position.x = mesh->mVertices[i].x;
        vertex.position.y = mesh->mVertices[i].y;
        vertex.position.z = mesh->mVertices[i].z;

        // Normal
        if (mesh->HasNormals()) {
            vertex.normal.x = mesh->mNormals[i].x;
            vertex.normal.y = mesh->mNormals[i].y;
            vertex.normal.z = mesh->mNormals[i].z;
        }

        // Texture coordinates
        if (mesh->mTextureCoords[0]) {
            vertex.texCoords.x = mesh->mTextureCoords[0][i].x;
            vertex.texCoords.y = mesh->mTextureCoords[0][i].y;

            // Tangent and bitangent for normal mapping
            if (mesh->HasTangentsAndBitangents()) {
                vertex.tangent.x = mesh->mTangents[i].x;
                vertex.tangent.y = mesh->mTangents[i].y;
                vertex.tangent.z = mesh->mTangents[i].z;

                vertex.bitangent.x = mesh->mBitangents[i].x;
                vertex.bitangent.y = mesh->mBitangents[i].y;
                vertex.bitangent.z = mesh->mBitangents[i].z;
            }
        }
        else {
            vertex.texCoords = glm::vec2(0.0f, 0.0f);
        }

        // Default bone indices and weights to prevent uninitialized values
        vertex.boneIDs = glm::ivec4(-1, -1, -1, -1);
        vertex.weights = glm::vec4(0.0f);

        result.vertices.push_back(vertex);
    }

    // Process indices
    for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
        aiFace face = mesh->mFaces[i];
        for (unsigned int j = 0; j < face.mNumIndices; j++) {
            result.indices.push_back(face.mIndices[j]);
        }
    }

    // Process bones if present
    if (mesh->HasBones()) {
        for (unsigned int i = 0; i < mesh->mNumBones; i++) {
            aiBone* bone = mesh->mBones[i];
            std::string boneName = bone->mName.C_Str();
            int boneID = 0;

            // If bone isn't already in the bone list
            if (boneMap.find(boneName) == boneMap.end()) {
                boneID = bones.size();
                Bone newBone;
                newBone.name = boneName;
                newBone.id = boneID;

                // Convert from Assimp matrix to GLM
                aiMatrix4x4 offset = bone->mOffsetMatrix;
                newBone.offset = glm::mat4(
                    offset.a1, offset.b1, offset.c1, offset.d1,
                    offset.a2, offset.b2, offset.c2, offset.d2,
                    offset.a3, offset.b3, offset.c3, offset.d3,
                    offset.a4, offset.b4, offset.c4, offset.d4
                );

                bones.push_back(newBone);
                boneMap[boneName] = boneID;
            }
            else {
                boneID = boneMap[boneName];
            }

            // Add bone weights to vertices
            for (unsigned int j = 0; j < bone->mNumWeights; j++) {
                aiVertexWeight weight = bone->mWeights[j];
                unsigned int vertexID = weight.mVertexId;
                float weightValue = weight.mWeight;

                // Find the first available slot in the vertex
                for (int k = 0; k < 4; k++) {
                    if (result.vertices[vertexID].weights[k] == 0.0f) {
                        result.vertices[vertexID].boneIDs[k] = boneID;
                        result.vertices[vertexID].weights[k] = weightValue;
                        break;
                    }
                }
            }
        }

        // Normalize weights
        for (auto& vertex : result.vertices) {
            float sum = vertex.weights.x + vertex.weights.y + vertex.weights.z + vertex.weights.w;
            if (sum > 0.0f) {
                vertex.weights /= sum;
            }
        }
    }

    // Process material
    if (mesh->mMaterialIndex >= 0) {
        aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];

        // Process material properties
        result.material = processMaterial(material);

        // Load diffuse textures
        std::vector<Texture> diffuseMaps = loadMaterialTextures(material, aiTextureType_DIFFUSE, "texture_diffuse", scene);
        result.textures.insert(result.textures.end(), diffuseMaps.begin(), diffuseMaps.end());

        // Load specular textures
        std::vector<Texture> specularMaps = loadMaterialTextures(material, aiTextureType_SPECULAR, "texture_specular", scene);
        result.textures.insert(result.textures.end(), specularMaps.begin(), specularMaps.end());

        // Load normal maps
        std::vector<Texture> normalMaps = loadMaterialTextures(material, aiTextureType_NORMALS, "texture_normal", scene);
        result.textures.insert(result.textures.end(), normalMaps.begin(), normalMaps.end());

        // Load height maps
        std::vector<Texture> heightMaps = loadMaterialTextures(material, aiTextureType_HEIGHT, "texture_height", scene);
        result.textures.insert(result.textures.end(), heightMaps.begin(), heightMaps.end());

        // Load emissive maps
        std::vector<Texture> emissiveMaps = loadMaterialTextures(material, aiTextureType_EMISSIVE, "texture_emissive", scene);
        result.textures.insert(result.textures.end(), emissiveMaps.begin(), emissiveMaps.end());

        // Set hasTexture flag
        result.material.hasTexture = !result.textures.empty();
    }

    // Create OpenGL buffers
    glGenVertexArrays(1, &result.VAO);
    glGenBuffers(1, &result.VBO);
    glGenBuffers(1, &result.EBO);

    glBindVertexArray(result.VAO);

    glBindBuffer(GL_ARRAY_BUFFER, result.VBO);
    glBufferData(GL_ARRAY_BUFFER, result.vertices.size() * sizeof(Vertex), &result.vertices[0], GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, result.EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, result.indices.size() * sizeof(unsigned int), &result.indices[0], GL_STATIC_DRAW);

    // Vertex positions
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);

    // Vertex normals
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));

    // Vertex texture coords
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, texCoords));

    // Vertex tangent
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tangent));

    // Vertex bitangent
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, bitangent));

    // Vertex bone IDs
    glEnableVertexAttribArray(5);
    glVertexAttribIPointer(5, 4, GL_INT, sizeof(Vertex), (void*)offsetof(Vertex, boneIDs));

    // Vertex weights
    glEnableVertexAttribArray(6);
    glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, weights));

    glBindVertexArray(0);

    return result;
}

Model::Material Model::processMaterial(aiMaterial* material) {
    Material result;

    // Set default values
    result.ambient = glm::vec3(0.2f, 0.2f, 0.2f);
    result.diffuse = glm::vec3(0.8f, 0.8f, 0.8f);
    result.specular = glm::vec3(0.5f, 0.5f, 0.5f);
    result.emissive = glm::vec3(0.0f, 0.0f, 0.0f);
    result.shininess = 32.0f;
    result.opacity = 1.0f;
    result.hasTexture = false;

    // Try to load material properties
    aiColor3D color(0.0f, 0.0f, 0.0f);
    float value = 0.0f;

    // Ambient color
    if (material->Get(AI_MATKEY_COLOR_AMBIENT, color) == AI_SUCCESS) {
        result.ambient = glm::vec3(color.r, color.g, color.b);
    }

    // Diffuse color
    if (material->Get(AI_MATKEY_COLOR_DIFFUSE, color) == AI_SUCCESS) {
        result.diffuse = glm::vec3(color.r, color.g, color.b);
    }

    // Specular color
    if (material->Get(AI_MATKEY_COLOR_SPECULAR, color) == AI_SUCCESS) {
        result.specular = glm::vec3(color.r, color.g, color.b);
    }

    // Emissive color
    if (material->Get(AI_MATKEY_COLOR_EMISSIVE, color) == AI_SUCCESS) {
        result.emissive = glm::vec3(color.r, color.g, color.b);
    }

    // Shininess
    if (material->Get(AI_MATKEY_SHININESS, value) == AI_SUCCESS) {
        result.shininess = value;
    }

    // Opacity
    if (material->Get(AI_MATKEY_OPACITY, value) == AI_SUCCESS) {
        result.opacity = value;
    }

    return result;
}

std::vector<Model::Texture> Model::loadMaterialTextures(aiMaterial* mat, aiTextureType type,
    const std::string& typeName, const aiScene* scene) {
    std::vector<Texture> textures;

    for (unsigned int i = 0; i < mat->GetTextureCount(type); i++) {
        aiString str;
        mat->GetTexture(type, i, &str);

        // Check if texture was loaded before
        bool skip = false;
        for (unsigned int j = 0; j < loadedTextures.size(); j++) {
            if (loadedTextures.find(str.C_Str()) != loadedTextures.end()) {
                Texture texture;
                texture.id = loadedTextures[str.C_Str()];
                texture.type = typeName;
                texture.path = str.C_Str();
                texture.embedded = false;
                textures.push_back(texture);
                skip = true;
                break;
            }
        }

        if (!skip) {
            Texture texture;
            texture.type = typeName;
            texture.path = str.C_Str();

            // Check if the texture is embedded
            const aiTexture* embeddedTexture = scene->GetEmbeddedTexture(str.C_Str());
            if (embeddedTexture) {
                texture.id = loadEmbeddedTexture(embeddedTexture);
                texture.embedded = true;
                std::cout << "Loaded embedded texture: " << str.C_Str() << std::endl;
            }
            else {
                // Check if texture file exists
                std::string fullPath = directory + "/" + str.C_Str();
                texture.id = loadTextureFromFile(fullPath);
                texture.embedded = false;
                std::cout << "Loaded texture: " << fullPath << std::endl;
            }

            textures.push_back(texture);
            loadedTextures[str.C_Str()] = texture.id;
        }
    }

    return textures;
}

unsigned int Model::loadTextureFromFile(const std::string& path) {
    unsigned int textureID;
    glGenTextures(1, &textureID);

    int width, height, nrComponents;
    unsigned char* data = stbi_load(path.c_str(), &width, &height, &nrComponents, 0);

    if (data) {
        GLenum format;
        if (nrComponents == 1) {
            format = GL_RED;
        }
        else if (nrComponents == 3) {
            format = GL_RGB;
        }
        else if (nrComponents == 4) {
            format = GL_RGBA;
        }
        else {
            format = GL_RGB;
            std::cout << "Unknown image format with " << nrComponents << " components" << std::endl;
        }

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);

        // Set texture parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        stbi_image_free(data);
    }
    else {
        std::cout << "Texture failed to load at path: " << path << std::endl;
        stbi_image_free(data);
        return 0;
    }

    return textureID;
}

unsigned int Model::loadEmbeddedTexture(const aiTexture* embeddedTexture) {
    unsigned int textureID;
    glGenTextures(1, &textureID);

    int width, height, nrComponents;
    unsigned char* data;

    // Check if texture is compressed (mHeight == 0 means it's compressed)
    if (embeddedTexture->mHeight == 0) {
        // Compressed texture format, use stbi to load it
        data = stbi_load_from_memory(
            reinterpret_cast<const stbi_uc*>(embeddedTexture->pcData),
            embeddedTexture->mWidth,
            &width, &height, &nrComponents, 0);

        if (data) {
            GLenum format;
            if (nrComponents == 1) {
                format = GL_RED;
            }
            else if (nrComponents == 3) {
                format = GL_RGB;
            }
            else if (nrComponents == 4) {
                format = GL_RGBA;
            }
            else {
                format = GL_RGB;
                std::cout << "Unknown embedded texture format with " << nrComponents << " components" << std::endl;
            }

            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
            glGenerateMipmap(GL_TEXTURE_2D);

            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

            stbi_image_free(data);
        }
        else {
            std::cout << "Embedded texture failed to load" << std::endl;
            stbi_image_free(data);
            return 0;
        }
    }
    else {
        // Uncompressed format, raw RGBA data
        width = embeddedTexture->mWidth;
        height = embeddedTexture->mHeight;

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, embeddedTexture->pcData);
        glGenerateMipmap(GL_TEXTURE_2D);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    }

    return textureID;
}

void Model::loadAnimations(const aiScene* scene) {
    if (!scene->HasAnimations()) {
        return;
    }

    // Process each animation
    for (unsigned int i = 0; i < scene->mNumAnimations; i++) {
        aiAnimation* anim = scene->mAnimations[i];
        Animation animation;

        // Store animation properties
        animation.name = anim->mName.length > 0 ? anim->mName.C_Str() : "Animation_" + std::to_string(i);
        animation.duration = anim->mDuration;
        animation.ticksPerSecond = anim->mTicksPerSecond != 0 ? anim->mTicksPerSecond : 25.0f;
        animation.loop = true;  // Default to looping

        std::cout << "Found animation " << i << ": " << animation.name
            << " (Duration: " << animation.duration
            << ", TPS: " << animation.ticksPerSecond
            << ", Channels: " << anim->mNumChannels << ")" << std::endl;

        animations.push_back(animation);
    }
}

void Model::update(float deltaTime) {
    if (animations.empty() || !scene || !scene->HasAnimations()) {
        return;  // No animations to update
    }

    Animation& animation = animations[currentAnimation];

    // Update animation time
    animationTime += deltaTime * animation.ticksPerSecond;

    // Handle looping
    if (animation.loop) {
        while (animationTime >= animation.duration) {
            animationTime -= animation.duration;
        }
    }
    else if (animationTime > animation.duration) {
        animationTime = animation.duration;
    }

    // Update bone transforms
    updateBoneTransforms(animationTime, scene);
}

void Model::updateBoneTransforms(float time, const aiScene* scene) {
    if (!scene->HasAnimations() || currentAnimation >= scene->mNumAnimations) {
        return;
    }

    aiAnimation* animation = scene->mAnimations[currentAnimation];

    // Create a map to store node transformations
    std::map<std::string, aiMatrix4x4> nodeTransforms;

    // Helper function to find animation channel for a node
    auto findNodeAnim = [animation](const std::string& nodeName) -> const aiNodeAnim* {
        for (unsigned int i = 0; i < animation->mNumChannels; i++) {
            if (animation->mChannels[i]->mNodeName.C_Str() == nodeName) {
                return animation->mChannels[i];
            }
        }
        return nullptr;
        };

    // Helper function to interpolate position keys
    auto interpolatePosition = [](float time, const aiNodeAnim* nodeAnim) -> aiVector3D {
        if (nodeAnim->mNumPositionKeys == 1) {
            return nodeAnim->mPositionKeys[0].mValue;
        }

        unsigned int index = 0;
        for (unsigned int i = 0; i < nodeAnim->mNumPositionKeys - 1; i++) {
            if (time < nodeAnim->mPositionKeys[i + 1].mTime) {
                index = i;
                break;
            }
        }

        unsigned int nextIndex = index + 1;
        float t = (time - nodeAnim->mPositionKeys[index].mTime) /
            (nodeAnim->mPositionKeys[nextIndex].mTime - nodeAnim->mPositionKeys[index].mTime);

        const aiVector3D& start = nodeAnim->mPositionKeys[index].mValue;
        const aiVector3D& end = nodeAnim->mPositionKeys[nextIndex].mValue;

        return start + (end - start) * t;
        };

    // Helper function to interpolate rotation keys
    auto interpolateRotation = [](float time, const aiNodeAnim* nodeAnim) -> aiQuaternion {
        if (nodeAnim->mNumRotationKeys == 1) {
            return nodeAnim->mRotationKeys[0].mValue;
        }

        unsigned int index = 0;
        for (unsigned int i = 0; i < nodeAnim->mNumRotationKeys - 1; i++) {
            if (time < nodeAnim->mRotationKeys[i + 1].mTime) {
                index = i;
                break;
            }
        }

        unsigned int nextIndex = index + 1;
        float t = (time - nodeAnim->mRotationKeys[index].mTime) /
            (nodeAnim->mRotationKeys[nextIndex].mTime - nodeAnim->mRotationKeys[index].mTime);

        const aiQuaternion& start = nodeAnim->mRotationKeys[index].mValue;
        const aiQuaternion& end = nodeAnim->mRotationKeys[nextIndex].mValue;

        aiQuaternion result;
        aiQuaternion::Interpolate(result, start, end, t);
        return result.Normalize();
        };

    // Helper function to interpolate scaling keys
    auto interpolateScaling = [](float time, const aiNodeAnim* nodeAnim) -> aiVector3D {
        if (nodeAnim->mNumScalingKeys == 1) {
            return nodeAnim->mScalingKeys[0].mValue;
        }

        unsigned int index = 0;
        for (unsigned int i = 0; i < nodeAnim->mNumScalingKeys - 1; i++) {
            if (time < nodeAnim->mScalingKeys[i + 1].mTime) {
                index = i;
                break;
            }
        }

        unsigned int nextIndex = index + 1;
        float t = (time - nodeAnim->mScalingKeys[index].mTime) /
            (nodeAnim->mScalingKeys[nextIndex].mTime - nodeAnim->mScalingKeys[index].mTime);

        const aiVector3D& start = nodeAnim->mScalingKeys[index].mValue;
        const aiVector3D& end = nodeAnim->mScalingKeys[nextIndex].mValue;

        return start + (end - start) * t;
        };

    // Recursive function to calculate node transformations
    std::function<void(aiNode*, const aiMatrix4x4&)> calculateNodeTransform;
    calculateNodeTransform = [&](aiNode* node, const aiMatrix4x4& parentTransform) {
        std::string nodeName = node->mName.C_Str();
        aiMatrix4x4 nodeTransformation = node->mTransformation;

        // Find animation channel for this node
        const aiNodeAnim* nodeAnim = findNodeAnim(nodeName);
        if (nodeAnim) {
            // Interpolate transformation
            aiVector3D position = interpolatePosition(time, nodeAnim);
            aiQuaternion rotation = interpolateRotation(time, nodeAnim);
            aiVector3D scaling = interpolateScaling(time, nodeAnim);

            // Create transformation matrix
            aiMatrix4x4 positionMat;
            aiMatrix4x4::Translation(position, positionMat);

            aiMatrix4x4 rotationMat = aiMatrix4x4(rotation.GetMatrix());

            aiMatrix4x4 scalingMat;
            aiMatrix4x4::Scaling(scaling, scalingMat);

            // Combine transformations
            nodeTransformation = positionMat * rotationMat * scalingMat;
        }

        // Calculate global transformation
        aiMatrix4x4 globalTransform = parentTransform * nodeTransformation;
        nodeTransforms[nodeName] = globalTransform;

        // Process children
        for (unsigned int i = 0; i < node->mNumChildren; i++) {
            calculateNodeTransform(node->mChildren[i], globalTransform);
        }
        };

    // Start from the root node with identity matrix
    aiMatrix4x4 identity;
    calculateNodeTransform(scene->mRootNode, identity);

    // Update bone transforms
    for (const auto& bone : bones) {
        if (nodeTransforms.find(bone.name) != nodeTransforms.end()) {
            aiMatrix4x4 globalTransform = nodeTransforms[bone.name];

            // Convert from Assimp matrix to GLM
            glm::mat4 offset = bone.offset;
            glm::mat4 globalTransformGLM = glm::mat4(
                globalTransform.a1, globalTransform.b1, globalTransform.c1, globalTransform.d1,
                globalTransform.a2, globalTransform.b2, globalTransform.c2, globalTransform.d2,
                globalTransform.a3, globalTransform.b3, globalTransform.c3, globalTransform.d3,
                globalTransform.a4, globalTransform.b4, globalTransform.c4, globalTransform.d4
            );

            // Final bone transform
            boneTransforms[bone.id] = globalTransformGLM * offset;
        }
    }
}

void Model::setAnimation(int animIndex) {
    if (animIndex >= 0 && animIndex < animations.size()) {
        currentAnimation = animIndex;
        animationTime = 0.0f; // Reset animation time when switching
        std::cout << "Switched to animation " << animIndex;
        if (scene && scene->HasAnimations() && animIndex < scene->mNumAnimations) {
            std::cout << " (" << scene->mAnimations[animIndex]->mName.C_Str() << ")";
        }
        std::cout << std::endl;
    }
}

int Model::getAnimationCount() const {
    return animations.size();
}

std::string Model::getAnimationName(int index) const {
    if (index >= 0 && index < animations.size()) {
        return animations[index].name;
    }
    return "Unknown";
}

bool Model::hasTextures() const {
    for (const auto& mesh : meshes) {
        if (!mesh.textures.empty()) {
            return true;
        }
    }
    return false;
}

void Model::setMeshMaterial(int meshIndex, const glm::vec3& ambient, const glm::vec3& diffuse,
    const glm::vec3& specular, float shininess) {
    if (meshIndex >= 0 && meshIndex < meshes.size()) {
        meshes[meshIndex].material.ambient = ambient;
        meshes[meshIndex].material.diffuse = diffuse;
        meshes[meshIndex].material.specular = specular;
        meshes[meshIndex].material.shininess = shininess;
    }
}

void Model::draw(std::shared_ptr<Program> shader) {
    // Set model matrix based on position, rotation and scale
    glm::mat4 model = glm::mat4(1.0f);
    model = glm::translate(model, position);
    model = glm::rotate(model, rotationAngle, rotationAxis);
    model = glm::scale(model, scale);

    glUniformMatrix4fv(shader->getUniform("M"), 1, GL_FALSE, glm::value_ptr(model));

    // Set bone transforms if needed
    if (!boneTransforms.empty()) {
        for (unsigned int i = 0; i < boneTransforms.size(); i++) {
            std::string uniformName = "boneTransforms[" + std::to_string(i) + "]";
            glUniformMatrix4fv(shader->getUniform(uniformName), 1, GL_FALSE, glm::value_ptr(boneTransforms[i]));
        }
    }

    // Draw meshes
    for (unsigned int i = 0; i < meshes.size(); i++) {
        const auto& mesh = meshes[i];

        // Set material properties
        glUniform3fv(shader->getUniform("material.ambient"), 1, glm::value_ptr(mesh.material.ambient));
        glUniform3fv(shader->getUniform("material.diffuse"), 1, glm::value_ptr(mesh.material.diffuse));
        glUniform3fv(shader->getUniform("material.specular"), 1, glm::value_ptr(mesh.material.specular));
        glUniform3fv(shader->getUniform("material.emissive"), 1, glm::value_ptr(mesh.material.emissive));
        glUniform1f(shader->getUniform("material.shininess"), mesh.material.shininess);
        glUniform1f(shader->getUniform("material.opacity"), mesh.material.opacity);

        // Bind appropriate textures
        unsigned int diffuseNr = 1;
        unsigned int specularNr = 1;
        unsigned int normalNr = 1;
        unsigned int heightNr = 1;
        unsigned int emissiveNr = 1;

        glUniform1i(shader->getUniform("material.hasTexture"), !mesh.textures.empty());

        for (unsigned int j = 0; j < mesh.textures.size(); j++) {
            // Activate proper texture unit before binding
            glActiveTexture(GL_TEXTURE0 + j);

            // Retrieve texture number (e.g. diffuse_1, diffuse_2, etc.)
            std::string number;
            std::string name = mesh.textures[j].type;

            if (name == "texture_diffuse") {
                number = std::to_string(diffuseNr++);
            }
            else if (name == "texture_specular") {
                number = std::to_string(specularNr++);
            }
            else if (name == "texture_normal") {
                number = std::to_string(normalNr++);
            }
            else if (name == "texture_height") {
                number = std::to_string(heightNr++);
            }
            else if (name == "texture_emissive") {
                number = std::to_string(emissiveNr++);
            }

            // Set the sampler to the correct texture unit
            glUniform1i(shader->getUniform(name + number), j);

            // Bind the texture
            glBindTexture(GL_TEXTURE_2D, mesh.textures[j].id);
        }

        // Draw mesh
        glBindVertexArray(mesh.VAO);
        glDrawElements(GL_TRIANGLES, mesh.indices.size(), GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        // Reset active texture
        glActiveTexture(GL_TEXTURE0);
    }
}